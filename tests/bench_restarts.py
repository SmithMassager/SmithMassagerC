#!/usr/bin/env python3
"""
Estimate the per-attempt success probability of the Las Vegas algorithm.

Runs smith-massager multiple times for each prime and counts the number
of Las Vegas attempts (restarts) per run.

Usage:
    python3 tests/bench_restarts.py [--binary PATH] [--trials N] [--max-n N]
"""

import argparse
import csv
import os
import statistics
import subprocess

DEFAULT_BINARY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "bin", "smith-massager"
)

PRIMES = [
    11, 23, 53, 97, 101, 127, 151, 199, 251, 307,
    401, 503, 601, 701, 809, 907, 1009,
]


def main():
    parser = argparse.ArgumentParser(description="Restart probability estimation")
    parser.add_argument("--binary", default=DEFAULT_BINARY)
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--max-n", type=int, default=1009)
    parser.add_argument("--output", default=None,
                        help="Output CSV path (default: stdout)")
    args = parser.parse_args()

    primes = [p for p in PRIMES if p <= args.max_n]

    print(f"Binary: {args.binary}")
    print(f"Trials per prime: {args.trials}")
    print()

    results = []  # (n, trial, attempts)

    for n in primes:
        print(f"n={n}: ", end="", flush=True)
        for t in range(1, args.trials + 1):
            try:
                r = subprocess.run(
                    [args.binary, str(n)],
                    capture_output=True, text=True, timeout=600
                )
                attempts = r.stdout.count("largestInvFactor:")
                results.append((n, t, attempts))
                print(f"{attempts} ", end="", flush=True)
            except subprocess.TimeoutExpired:
                results.append((n, t, -1))
                print("T ", end="", flush=True)
        print()

    # Write CSV
    if args.output:
        with open(args.output, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["n", "trial", "attempts"])
            for row in results:
                w.writerow(row)
        print(f"\nSaved to {args.output}")

    # Print summary
    print(f"\n{'n':>6} {'trials':>6} {'mean_att':>10} {'med_att':>8} "
          f"{'min':>5} {'max':>5} {'P(success)':>12}")
    print("-" * 60)

    from collections import defaultdict
    by_n = defaultdict(list)
    for n, t, att in results:
        if att > 0:
            by_n[n].append(att)

    for n in primes:
        if not by_n[n]:
            print(f"{n:>6} {'---':>6}")
            continue
        atts = by_n[n]
        mean_a = statistics.mean(atts)
        med_a = statistics.median(atts)
        # P(success) = total_successes / total_attempts
        # Each run has exactly 1 success and (attempts-1) failures
        total_att = sum(atts)
        total_succ = len(atts)
        p_success = total_succ / total_att
        print(f"{n:>6} {len(atts):>6} {mean_a:>10.2f} {med_a:>8.0f} "
              f"{min(atts):>5} {max(atts):>5} {p_success:>12.3f}")


if __name__ == "__main__":
    main()
