"""Specialized banded Gaussian elimination for SPAD solver linear systems."""

from __future__ import annotations


def row_echelon(matrix: list[float], bins: int) -> None:
    columns = bins

    for i in range(columns - 1):
        mult_factor = matrix[(i + 1) * (columns + 1) + i] / matrix[i * (columns + 1) + i]

        matrix[(i + 1) * (columns + 1) + i] = 0.0
        matrix[(i + 1) * (columns + 1) + bins - 1] -= (
            mult_factor * matrix[i * (columns + 1) + bins - 1]
        )
        matrix[(i + 1) * (columns + 1) + bins] -= mult_factor * matrix[i * (columns + 1) + bins]

    for i in range(columns - 1):
        mult_factor = matrix[i * (columns + 1) + bins - 1] / matrix[(bins - 1) * (columns + 1) + bins - 1]

        matrix[i * (columns + 1) + bins] -= mult_factor * matrix[(bins - 1) * (columns + 1) + bins]
        matrix[i * (columns + 1) + bins - 1] = 0.0

    for i in range(columns):
        matrix[i * (columns + 1) + bins] /= matrix[i * (columns + 1) + i]
        matrix[i * (columns + 1) + i] = 1.0


def sol_vector(matrix: list[float], bins: int) -> list[float]:
    return [matrix[i * (bins + 1) + bins] for i in range(bins)]


def solve_system(matrix: list[float], bins: int) -> list[float]:
    row_echelon(matrix, bins)
    return sol_vector(matrix, bins)
