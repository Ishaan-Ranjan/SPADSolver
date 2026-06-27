"""Tests for interpolation functions using shared test vectors."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from spadsolver import alpha, beta, gamma

_SHARED = Path(__file__).resolve().parents[2] / "shared" / "test_vectors"

K_FIELD = 3.0e6
K_TEMP = 298.0


def _load(name: str) -> dict:
    with (_SHARED / name).open(encoding="utf-8") as handle:
        return json.load(handle)


@pytest.mark.parametrize(
    "case",
    _load("alpha_beta.json")["cases"],
    ids=lambda c: c["name"],
)
def test_alpha_beta_shared_vectors(case: dict) -> None:
    defaults = _load("alpha_beta.json")["defaults"]
    fn = alpha if case["fn"] == "alpha" else beta
    actual = fn(case["al_frac"], defaults["field"], defaults["temp"])
    assert abs(actual - case["expected"]) <= defaults["tolerance"]


@pytest.mark.parametrize(
    "case",
    _load("gamma.json")["cases"],
    ids=lambda c: c["name"],
)
def test_gamma_shared_vectors(case: dict) -> None:
    meta = _load("gamma.json")
    actual = gamma(case["al_frac"], case["energy"])
    assert abs(actual - case["expected"]) <= meta["tolerance"]


def test_alpha_decreases_with_al_fraction() -> None:
    assert alpha(0.1, K_FIELD, K_TEMP) > alpha(0.5, K_FIELD, K_TEMP)


def test_beta_is_positive_for_typical_inputs() -> None:
    value = beta(0.3, K_FIELD, K_TEMP)
    assert value > 0.0
    assert math.isfinite(value)


def test_gamma_is_positive() -> None:
    value = gamma(0.38, 3.75)
    assert value > 0.0
    assert math.isfinite(value)
