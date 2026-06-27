#!/usr/bin/env python3
"""Benchmark serial vs Ray performance for parameter sweeps and linear solves.

This is a standalone harness (not run by pytest). Install Ray first:

    pip install -r requirements-ray.txt
    python benchmarks/benchmark_ray.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import statistics
import sys
import time
from dataclasses import asdict, dataclass
from typing import Any

try:
    import ray
except ImportError:
    print("Ray is not installed. Run: pip install -r requirements-ray.txt", file=sys.stderr)
    sys.exit(1)

from spadsolver.ray_tasks import (
    evaluate_batch,
    parallel_parameter_sweep,
    parallel_solve_systems,
    solve_linear_system,
)

ParameterPoint = tuple[float, float, float, float]


@dataclass
class TimingResult:
    label: str
    seconds: float
    num_items: int

    @property
    def items_per_second(self) -> float:
        if self.seconds <= 0:
            return float("inf")
        return self.num_items / self.seconds


@dataclass
class BenchmarkReport:
    parameter_sweep: dict[str, Any]
    linear_solve: dict[str, Any]


def build_parameter_grid(num_al_fracs: int, num_energies: int) -> list[ParameterPoint]:
    al_fracs = [i / max(num_al_fracs - 1, 1) * 0.5 for i in range(num_al_fracs)]
    energies = [2.5 + i * (4.0 - 2.5) / max(num_energies - 1, 1) for i in range(num_energies)]
    field, temp = 3.0e6, 298.0
    return [(a, field, temp, e) for a, e in itertools.product(al_fracs, energies)]


def build_linear_systems(count: int) -> list[tuple[list[list[float]], list[float]]]:
    systems = [
        ([[2.0, 1.0], [1.0, 3.0]], [5.0, 9.0]),
        ([[4.0, 0.0, 1.0], [2.0, 2.0, 0.0], [0.0, 3.0, 1.0]], [9.0, 6.0, 7.0]),
        ([[5.0]], [10.0]),
    ]
    return [systems[i % len(systems)] for i in range(count)]


def time_call(label: str, num_items: int, func, *args, **kwargs) -> TimingResult:
    start = time.perf_counter()
    func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return TimingResult(label=label, seconds=elapsed, num_items=num_items)


def run_timed(repeats: int, warmup: int, label: str, num_items: int, func, *args, **kwargs) -> TimingResult:
    for _ in range(warmup):
        func(*args, **kwargs)

    samples = []
    for _ in range(repeats):
        result = time_call(label, num_items, func, *args, **kwargs)
        samples.append(result.seconds)

    median_seconds = statistics.median(samples)
    return TimingResult(label=label, seconds=median_seconds, num_items=num_items)


def benchmark_parameter_sweep(
    points: list[ParameterPoint],
    *,
    chunk_size: int,
    num_cpus: int | None,
    repeats: int,
    warmup: int,
) -> dict[str, Any]:
    serial = run_timed(repeats, warmup, "serial", len(points), evaluate_batch, points)
    parallel = run_timed(
        repeats,
        warmup,
        "ray",
        len(points),
        parallel_parameter_sweep,
        points,
        chunk_size=chunk_size,
        num_cpus=num_cpus,
    )
    speedup = serial.seconds / parallel.seconds if parallel.seconds > 0 else float("inf")
    return {
        "num_points": len(points),
        "chunk_size": chunk_size,
        "num_cpus": num_cpus,
        "serial_seconds": serial.seconds,
        "ray_seconds": parallel.seconds,
        "speedup": speedup,
        "serial_items_per_second": serial.items_per_second,
        "ray_items_per_second": parallel.items_per_second,
    }


def benchmark_linear_solve(
    systems: list[tuple[list[list[float]], list[float]]],
    *,
    num_cpus: int | None,
    repeats: int,
    warmup: int,
) -> dict[str, Any]:
    def serial_solve() -> None:
        for coeffs, rhs in systems:
            solve_linear_system(coeffs, rhs)

    serial = run_timed(repeats, warmup, "serial", len(systems), serial_solve)
    parallel = run_timed(
        repeats,
        warmup,
        "ray",
        len(systems),
        parallel_solve_systems,
        systems,
        num_cpus=num_cpus,
    )
    speedup = serial.seconds / parallel.seconds if parallel.seconds > 0 else float("inf")
    return {
        "num_systems": len(systems),
        "num_cpus": num_cpus,
        "serial_seconds": serial.seconds,
        "ray_seconds": parallel.seconds,
        "speedup": speedup,
        "serial_items_per_second": serial.items_per_second,
        "ray_items_per_second": parallel.items_per_second,
    }


def print_report(report: BenchmarkReport) -> None:
    sweep = report.parameter_sweep
    solve = report.linear_solve

    print("=" * 60)
    print("SPADSolver Ray performance benchmark")
    print("=" * 60)
    print()
    print("Parameter sweep (alpha / beta / gamma)")
    print(f"  Points:              {sweep['num_points']}")
    print(f"  Ray chunk_size:      {sweep['chunk_size']}")
    print(f"  Ray num_cpus:        {sweep['num_cpus']}")
    print(f"  Serial (median):     {sweep['serial_seconds']:.4f}s  "
          f"({sweep['serial_items_per_second']:.1f} points/s)")
    print(f"  Ray (median):        {sweep['ray_seconds']:.4f}s  "
          f"({sweep['ray_items_per_second']:.1f} points/s)")
    print(f"  Speedup (serial/ray): {sweep['speedup']:.2f}x")
    print()
    print("Linear solve batch")
    print(f"  Systems:             {solve['num_systems']}")
    print(f"  Ray num_cpus:        {solve['num_cpus']}")
    print(f"  Serial (median):     {solve['serial_seconds']:.4f}s  "
          f"({solve['serial_items_per_second']:.1f} systems/s)")
    print(f"  Ray (median):        {solve['ray_seconds']:.4f}s  "
          f"({solve['ray_items_per_second']:.1f} systems/s)")
    print(f"  Speedup (serial/ray): {solve['speedup']:.2f}x")
    print()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--al-fracs", type=int, default=20, help="Grid points along Al fraction")
    parser.add_argument("--energies", type=int, default=20, help="Grid points along energy")
    parser.add_argument("--systems", type=int, default=200, help="Number of linear systems to solve")
    parser.add_argument("--chunk-size", type=int, default=64, help="Points per Ray task")
    parser.add_argument("--cpus", type=int, default=None, help="Ray worker CPU limit")
    parser.add_argument("--repeats", type=int, default=3, help="Timed repetitions (median reported)")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs before timing")
    parser.add_argument("--json", action="store_true", help="Print JSON instead of text table")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    points = build_parameter_grid(args.al_fracs, args.energies)
    systems = build_linear_systems(args.systems)

    ray.init(num_cpus=args.cpus, ignore_reinit_error=True)

    report = BenchmarkReport(
        parameter_sweep=benchmark_parameter_sweep(
            points,
            chunk_size=args.chunk_size,
            num_cpus=args.cpus,
            repeats=args.repeats,
            warmup=args.warmup,
        ),
        linear_solve=benchmark_linear_solve(
            systems,
            num_cpus=args.cpus,
            repeats=args.repeats,
            warmup=args.warmup,
        ),
    )

    if args.json:
        print(json.dumps(asdict(report), indent=2))
    else:
        print_report(report)


if __name__ == "__main__":
    main()
