"""Ionization and absorption coefficient interpolation."""

from __future__ import annotations

import functools
import json
import math
from pathlib import Path

_SHARED_CONSTANTS = (
    Path(__file__).resolve().parents[2] / "shared" / "constants" / "interpolation_tables.json"
)


@functools.lru_cache(maxsize=1)
def _load_tables() -> dict:
    with _SHARED_CONSTANTS.open(encoding="utf-8") as handle:
        return json.load(handle)


def load_interpolation_tables() -> dict:
    """Return a copy of the shared interpolation lookup tables."""
    return dict(_load_tables())


def lininterpolate(x: float, x1: float, y1: float, x2: float, y2: float) -> float:
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)


def _gan_ionization(front: float, top: float, field: float) -> float:
    return front * math.exp(-top / field)


def _front_gan(a: float, b: float, temp: float) -> float:
    return a * (1.0 + b * (temp - 298.0))


def _top_gan(c: float, d: float, temp: float) -> float:
    return -1.0 * c * (1.0 + d * (temp - 298.0))


def _algan_ionization(a: float, b: float, field: float) -> float:
    return a * math.exp(-b / field)


def alpha(al_frac: float, field: float, temp: float) -> float:
    a_n = 2.69e7
    b_n = 2.00e-3
    c_n = 2.27e7
    d_n = 5.00e-4

    front_gan_n = _front_gan(a_n, b_n, temp)
    top_gan_n = _top_gan(c_n, d_n, temp)

    front_alpha = 7.82e6
    top_alpha = 3.7e7

    return lininterpolate(
        al_frac,
        0.0,
        _gan_ionization(front_gan_n, top_gan_n, field),
        0.65,
        _algan_ionization(front_alpha, top_alpha, field),
    )


def beta(al_frac: float, field: float, temp: float) -> float:
    a_p = 4.32e6
    b_p = 2.00e-3
    c_p = 1.31e7
    d_p = 9.00e-4

    front_gan_p = _front_gan(a_p, b_p, temp)
    top_gan_p = _top_gan(c_p, d_p, temp)

    front_beta = 5.65e4
    top_beta = 7.04e6

    return lininterpolate(
        al_frac,
        0.0,
        _gan_ionization(front_gan_p, top_gan_p, field),
        0.65,
        _algan_ionization(front_beta, top_beta, field),
    )


def _gamma_exp_at_index(tables: dict, index: int, energy: float) -> float:
    sampled_energies = tables["sampled_energies"][index]
    sampled_abs_exps = tables["sampled_abs_exps"][index]
    sampled_abs_asy = tables["sampled_abs_asy"][index]

    gamma_exp = 0.0
    for j in range(3):
        if energy <= sampled_energies[0]:
            gamma_exp = lininterpolate(
                energy,
                sampled_energies[j],
                sampled_abs_exps[j],
                sampled_energies[j + 1],
                sampled_abs_exps[j + 1],
            )
            break
        if sampled_energies[j] <= energy <= sampled_energies[j + 1]:
            gamma_exp = lininterpolate(
                energy,
                sampled_energies[j],
                sampled_abs_exps[j],
                sampled_energies[j + 1],
                sampled_abs_exps[j + 1],
            )
            break
        gamma_exp = sampled_abs_asy
    return gamma_exp


def gamma(al_frac: float, energy: float) -> float:
    tables = _load_tables()
    mole_fracs = tables["mole_fracs"]

    gamma_final = 0.0
    for i in range(5):
        if mole_fracs[i] <= al_frac <= mole_fracs[i + 1]:
            mol_low = mole_fracs[i]
            mol_high = mole_fracs[i + 1]
        elif i == 4:
            mol_low = mole_fracs[4]
            mol_high = mole_fracs[5]
        else:
            continue

        gamma_exp1 = _gamma_exp_at_index(tables, i, energy)
        gamma_exp2 = _gamma_exp_at_index(tables, i + 1, energy)
        gamma_final = lininterpolate(al_frac, mol_low, gamma_exp1, mol_high, gamma_exp2)
        gamma_final = 10.0**gamma_final
        break
    return gamma_final
