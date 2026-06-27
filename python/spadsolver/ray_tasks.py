"""Optional Ray wrappers for parallel parameter sweeps and batched solves."""

from __future__ import annotations

from typing import Any

try:
    import ray
except ImportError:
    ray = None  # type: ignore[assignment,misc]

from spadsolver import alpha, beta, gamma
from spadsolver.gaussian_elimination import solve_system

# (al_frac, field_V_per_cm, temp_K, energy_eV)
ParameterPoint = tuple[float, float, float, float]


def _require_ray() -> Any:
    if ray is None:
        raise ImportError(
            "Ray is not installed. Install with: pip install -e '.[ray]' from the python/ directory."
        )
    return ray


def evaluate_point(al_frac: float, field: float, temp: float, energy: float) -> dict[str, float]:
    """Evaluate alpha, beta, and gamma at a single parameter point."""
    return {
        "al_frac": al_frac,
        "field": field,
        "temp": temp,
        "energy": energy,
        "alpha": alpha(al_frac, field, temp),
        "beta": beta(al_frac, field, temp),
        "gamma": gamma(al_frac, energy),
    }


def evaluate_batch(points: list[ParameterPoint]) -> list[dict[str, float]]:
    """Evaluate many parameter points serially."""
    return [evaluate_point(al_frac, field, temp, energy) for al_frac, field, temp, energy in points]


def _make_augmented_matrix(coeffs: list[list[float]], rhs: list[float]) -> list[float]:
    bins = len(rhs)
    matrix = [0.0] * (bins * (bins + 1))
    for i in range(bins):
        for j in range(bins):
            matrix[i * (bins + 1) + j] = coeffs[i][j]
        matrix[i * (bins + 1) + bins] = rhs[i]
    return matrix


def solve_linear_system(coeffs: list[list[float]], rhs: list[float]) -> list[float]:
    """Solve a linear system using the banded Gaussian elimination routine."""
    bins = len(rhs)
    matrix = _make_augmented_matrix(coeffs, rhs)
    return solve_system(matrix, bins)


def parallel_parameter_sweep(
    points: list[ParameterPoint],
    *,
    chunk_size: int = 64,
    num_cpus: int | None = None,
) -> list[dict[str, float]]:
    """Evaluate parameter points in parallel using Ray batched tasks."""
    _require_ray()
    if not points:
        return []

    if num_cpus is not None:
        ray.init(num_cpus=num_cpus, ignore_reinit_error=True)
    else:
        ray.init(ignore_reinit_error=True)

    remote_batch = ray.remote(evaluate_batch)
    chunks = [points[i : i + chunk_size] for i in range(0, len(points), chunk_size)]
    futures = [remote_batch.remote(chunk) for chunk in chunks]
    nested = ray.get(futures)
    return [result for chunk in nested for result in chunk]


def parallel_solve_systems(
    systems: list[tuple[list[list[float]], list[float]]],
    *,
    num_cpus: int | None = None,
) -> list[list[float]]:
    """Solve many independent linear systems in parallel using Ray."""
    _require_ray()
    if not systems:
        return []

    if num_cpus is not None:
        ray.init(num_cpus=num_cpus, ignore_reinit_error=True)
    else:
        ray.init(ignore_reinit_error=True)

    remote_solve = ray.remote(solve_linear_system)
    futures = [remote_solve.remote(coeffs, rhs) for coeffs, rhs in systems]
    return ray.get(futures)


if ray is not None:

    @ray.remote
    def evaluate_point_remote(
        al_frac: float, field: float, temp: float, energy: float
    ) -> dict[str, float]:
        return evaluate_point(al_frac, field, temp, energy)

    @ray.remote
    def evaluate_batch_remote(points: list[ParameterPoint]) -> list[dict[str, float]]:
        return evaluate_batch(points)

    @ray.remote
    def solve_linear_system_remote(
        coeffs: list[list[float]], rhs: list[float]
    ) -> list[float]:
        return solve_linear_system(coeffs, rhs)

    @ray.remote
    class InterpolationActor:
        """Ray actor that loads interpolation tables once per worker."""

        def __init__(self) -> None:
            from spadsolver.interpolation import _load_tables

            self._tables = _load_tables()

        def gamma(self, al_frac: float, energy: float) -> float:
            from spadsolver.interpolation import gamma

            return gamma(al_frac, energy)

        def evaluate_point(
            self, al_frac: float, field: float, temp: float, energy: float
        ) -> dict[str, float]:
            return evaluate_point(al_frac, field, temp, energy)

else:
    evaluate_point_remote = None
    evaluate_batch_remote = None
    solve_linear_system_remote = None
    InterpolationActor = None
