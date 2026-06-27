"""Tests for Gaussian elimination using shared test vectors."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from spadsolver.gaussian_elimination import row_echelon, sol_vector, solve_system

_SHARED = Path(__file__).resolve().parents[2] / "shared" / "test_vectors"


def _load(name: str) -> dict:
    with (_SHARED / name).open(encoding="utf-8") as handle:
        return json.load(handle)


def _make_augmented_matrix(coeffs: list[list[float]], rhs: list[float]) -> list[float]:
    bins = len(rhs)
    matrix = [0.0] * (bins * (bins + 1))
    for i in range(bins):
        for j in range(bins):
            matrix[i * (bins + 1) + j] = coeffs[i][j]
        matrix[i * (bins + 1) + bins] = rhs[i]
    return matrix


def _residual(coeffs: list[list[float]], rhs: list[float], sol: list[float]) -> float:
    max_residual = 0.0
    for i, row in enumerate(coeffs):
        residual = -rhs[i]
        for j, coeff in enumerate(row):
            residual += coeff * sol[j]
        max_residual = max(max_residual, abs(residual))
    return max_residual


@pytest.mark.parametrize(
    "case",
    _load("gaussian_elimination.json")["cases"],
    ids=lambda c: c["name"],
)
def test_solve_shared_vectors(case: dict) -> None:
    meta = _load("gaussian_elimination.json")
    bins = len(case["rhs"])
    matrix = _make_augmented_matrix(case["coeffs"], case["rhs"])
    sol = solve_system(matrix, bins)

    for i, expected in enumerate(case["expected"]):
        assert abs(sol[i] - expected) <= meta["tolerance"]
    assert _residual(case["coeffs"], case["rhs"], sol) <= meta["tolerance"]


def test_sol_vector_extracts_last_column() -> None:
    matrix = [1.0, 0.0, 3.0, 0.0, 1.0, 5.0]
    sol = sol_vector(matrix, 2)
    assert sol == pytest.approx([3.0, 5.0])


def test_row_echelon_reduced_form_2x2() -> None:
    matrix = [2.0, 1.0, 5.0, 1.0, 3.0, 9.0]
    row_echelon(matrix, 2)
    assert matrix[0] == pytest.approx(1.0)
    assert matrix[1] == pytest.approx(0.0)
    assert matrix[2] == pytest.approx(1.2)
    assert matrix[3] == pytest.approx(0.0)
    assert matrix[4] == pytest.approx(1.0)
    assert matrix[5] == pytest.approx(2.6)
