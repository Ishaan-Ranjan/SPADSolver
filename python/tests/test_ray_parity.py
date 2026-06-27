"""Correctness tests: Ray parallel paths must match serial implementations."""

from __future__ import annotations

import itertools

import pytest

pytest.importorskip("ray")

from spadsolver.ray_tasks import (
    evaluate_batch,
    parallel_parameter_sweep,
    parallel_solve_systems,
    solve_linear_system,
)

pytestmark = pytest.mark.ray


def _assert_parameter_results_match(
    serial: list[dict[str, float]],
    parallel: list[dict[str, float]],
) -> None:
    assert len(serial) == len(parallel)
    keys = ("al_frac", "field", "temp", "energy", "alpha", "beta", "gamma")
    for s, p in zip(serial, parallel):
        for key in keys:
            assert s[key] == pytest.approx(p[key])


def test_parameter_sweep_matches_serial_small_grid() -> None:
    points = [
        (0.0, 3.0e6, 298.0, 3.0),
        (0.1, 3.0e6, 298.0, 3.125),
        (0.2, 3.0e6, 298.0, 3.5),
        (0.5, 3.0e6, 298.0, 4.0),
    ]
    serial = evaluate_batch(points)
    parallel = parallel_parameter_sweep(points, chunk_size=2, num_cpus=2)
    _assert_parameter_results_match(serial, parallel)


def test_parameter_sweep_matches_serial_product_grid() -> None:
    al_fracs = [0.0, 0.11, 0.20, 0.38]
    energies = [2.5, 3.0, 3.5, 4.0]
    field, temp = 3.0e6, 298.0
    points = [(a, field, temp, e) for a, e in itertools.product(al_fracs, energies)]

    serial = evaluate_batch(points)
    parallel = parallel_parameter_sweep(points, chunk_size=4, num_cpus=2)
    _assert_parameter_results_match(serial, parallel)


def test_parameter_sweep_matches_serial_single_chunk() -> None:
    """All work in one Ray task."""
    points = [(0.3, 3.0e6, 298.0, 3.6), (0.4, 3.0e6, 310.0, 3.8)]
    serial = evaluate_batch(points)
    parallel = parallel_parameter_sweep(points, chunk_size=64, num_cpus=2)
    _assert_parameter_results_match(serial, parallel)


def test_parameter_sweep_matches_serial_one_task_per_point() -> None:
    """Each point in its own Ray task."""
    points = [(al, 3.0e6, 298.0, 3.0) for al in (0.0, 0.05, 0.15, 0.25)]
    serial = evaluate_batch(points)
    parallel = parallel_parameter_sweep(points, chunk_size=1, num_cpus=2)
    _assert_parameter_results_match(serial, parallel)


def test_solve_systems_match_serial() -> None:
    systems = [
        ([[5.0]], [10.0]),
        ([[2.0, 1.0], [1.0, 3.0]], [5.0, 9.0]),
        ([[4.0, 0.0, 1.0], [2.0, 2.0, 0.0], [0.0, 3.0, 1.0]], [9.0, 6.0, 7.0]),
    ]
    serial = [solve_linear_system(coeffs, rhs) for coeffs, rhs in systems]
    parallel = parallel_solve_systems(systems, num_cpus=2)

    assert len(serial) == len(parallel)
    for s, p in zip(serial, parallel):
        assert s == pytest.approx(p)
