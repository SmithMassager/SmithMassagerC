#!/usr/bin/env python3
"""
Timing benchmarks for the Smith massager.

Runs smith-massager for a range of primes, records wall-clock time
and number of Las Vegas attempts per trial.

Usage:
    python3 tests/bench_timing.py [--binary PATH] [--trials N] [--max-n N]
"""

import argparse
import csv
import os
import statistics
import subprocess
import time

DEFAULT_BINARY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "bin", "smith-massager"
)

PRIMES = [
    53, 97, 101, 127, 151, 199, 251, 307, 401, 503,
    601, 701, 809, 907, 1009, 1499, 2027,
]


def run_trial(binary, n, timeout=600):
    """Run one trial, return (seconds, attempts)."""
    t0 = time.perf_counter()
    r = subprocess.run(
        [binary, str(n)],
        capture_output=True, text=True, timeout=timeout
    )
    elapsed = time.perf_counter() - t0

    if r.returncode != 0:
        return elapsed, -1

    attempts = r.stdout.count("largestInvFactor:")
    return elapsed, attempts


def main():
    parser = argparse.ArgumentParser(description="Timing benchmarks")
    parser.add_argument("--binary", default=DEFAULT_BINARY)
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=2027)
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: stdout)")
    args = parser.parse_args()

    primes = [p for p in PRIMES if p <= args.max_n]

    print(f"Binary: {args.binary}")
    print(f"Trials: {args.trials}")
    print(f"Primes: {primes}")
    print()

    results = []  # (n, trial, seconds, attempts)

    for n in primes:
        print(f"n={n}: ", end="", flush=True)
        for t in range(1, args.trials + 1):
            try:
                elapsed, attempts = run_trial(args.binary, n)
                results.append((n, t, elapsed, attempts))
                print(f"{elapsed:.2f}s({attempts}) ", end="", flush=True)
            except subprocess.TimeoutExpired:
                results.append((n, t, -1, -1))
                print("TIMEOUT ", end="", flush=True)
        print()

    # Write CSV
    if args.output:
        with open(args.output, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["n", "trial", "seconds", "attempts"])
            for row in results:
                w.writerow(row)
        print(f"\nSaved to {args.output}")

    # Print summary
    print(f"\n{'n':>6} {'trials':>6} {'med_time':>10} {'med_att':>8} "
          f"{'med_1pass':>10} {'min_time':>10} {'max_time':>10}")
    print("-" * 70)

    from collections import defaultdict
    by_n = defaultdict(list)
    for n, t, sec, att in results:
        if sec > 0 and att > 0:
            by_n[n].append((sec, att))

    for n in primes:
        if not by_n[n]:
            print(f"{n:>6} {'---':>6}")
            continue
        times = [s for s, a in by_n[n]]
        atts = [a for s, a in by_n[n]]
        one_pass = [s / a for s, a in by_n[n]]
        med_t = statistics.median(times)
        med_a = statistics.median(atts)
        med_1p = statistics.median(one_pass)
        print(f"{n:>6} {len(times):>6} {med_t:>10.3f} {med_a:>8.0f} "
              f"{med_1p:>10.3f} {min(times):>10.3f} {max(times):>10.3f}")


if __name__ == "__main__":
    main()
