"""Tests for optional Ray integration."""

from __future__ import annotations

import pytest

from spadsolver.ray_tasks import (
    evaluate_batch,
    evaluate_point,
    parallel_parameter_sweep,
    solve_linear_system,
)


def test_evaluate_point_serial() -> None:
    result = evaluate_point(0.0, 3.0e6, 298.0, 3.125)
    assert result["alpha"] > 0.0
    assert result["beta"] > 0.0
    assert result["gamma"] > 0.0


def test_evaluate_batch_serial() -> None:
    points = [(0.0, 3.0e6, 298.0, 3.0), (0.2, 3.0e6, 298.0, 3.5)]
    results = evaluate_batch(points)
    assert len(results) == 2


def test_solve_linear_system() -> None:
    sol = solve_linear_system([[5.0]], [10.0])
    assert sol == pytest.approx([2.0])


@pytest.mark.ray
def test_parallel_parameter_sweep() -> None:
    pytest.importorskip("ray")
    points = [(al, 3.0e6, 298.0, 3.0) for al in (0.0, 0.1, 0.2, 0.3)]
    results = parallel_parameter_sweep(points, chunk_size=2, num_cpus=2)
    assert len(results) == 4
    assert all("gamma" in row for row in results)
