#!/usr/bin/env python3
"""Example: parallel parameter sweep over alpha, beta, and gamma using Ray."""

from __future__ import annotations

import argparse
import itertools
import time

from spadsolver.ray_tasks import parallel_parameter_sweep


def build_grid(
    al_fracs: list[float],
    energies: list[float],
    field: float,
    temp: float,
) -> list[tuple[float, float, float, float]]:
    return [(al_frac, field, temp, energy) for al_frac, energy in itertools.product(al_fracs, energies)]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpus", type=int, default=None, help="Ray worker CPU limit")
    parser.add_argument("--chunk-size", type=int, default=64, help="Points per Ray task")
    args = parser.parse_args()

    al_fracs = [0.0, 0.11, 0.20, 0.38, 0.5]
    energies = [2.5, 3.0, 3.5, 4.0]
    field = 3.0e6
    temp = 298.0

    points = build_grid(al_fracs, energies, field, temp)
    print(f"Sweeping {len(points)} parameter points with Ray...")

    start = time.perf_counter()
    results = parallel_parameter_sweep(points, chunk_size=args.chunk_size, num_cpus=args.cpus)
    elapsed = time.perf_counter() - start

    print(f"Completed in {elapsed:.3f}s")
    print("First result:", results[0])
    print("Last result:", results[-1])


if __name__ == "__main__":
    main()
