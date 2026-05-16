#!/usr/bin/env python3
"""
Verify SmithMassagerC output against SymPy.

Usage:
    python3 verify_sympy.py [max_n] [--binary PATH]

Constructs the test matrix A[i,j] = i^j mod n for each prime n,
runs bin/smith-massager to get (U, M, T, S), and checks:

  1. S is in Smith form (s_1 | s_2 | ... | s_n)
  2. S matches SymPy's smith_normal_form
  3. B is integral, where B = diag(A, I) * [I; U, I] * [I, M; 0, T] * [I; 0, S^{-1}]
  4. |det B| = 1 (unimodularity)

Requires: sympy  (pip install sympy)
"""

import argparse
import os
import re
import subprocess
import sys
from fractions import Fraction

from sympy import Matrix, ZZ, eye, zeros, det
from sympy.matrices.normalforms import smith_normal_form as snf


DEFAULT_BINARY = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "bin", "smith-massager"
)

PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]


def build_test_matrix(n):
    """Build the n x n matrix A[i,j] = i^j mod n."""
    return Matrix(n, n, lambda i, j: pow(i, j, n))


def parse_output(stdout, n):
    """Parse smith-massager stdout into (U, M, T, S, success).

    The output format is:
        <n> <n>
        <timing lines...>
        U, M, T, S:
        [[...] [...] ...]   <- U (n x n)
        [[...] [...] ...]   <- M (n x n)
        [[...] [...] ...]   <- T (n x n)
        [[...] [...] ...]   <- S (n x 1)
        success: 1
    """
    lines = stdout.strip().split("\n")
    success = any(line.strip() == "success: 1" for line in lines)

    # Find the "U, M, T, S:" marker
    marker = -1
    for i, line in enumerate(lines):
        if line.strip().startswith("U, M, T, S"):
            marker = i
            break
    if marker == -1:
        return None, None, None, None, success

    # Collect matrix text after marker
    matrix_text = "\n".join(lines[marker + 1:])

    # Split into four matrices by finding top-level [[ ]]
    # Each matrix starts with [[ and ends with ]
    matrices = []
    current = []
    depth = 0
    for char in matrix_text:
        if char == '[':
            depth += 1
        elif char == ']':
            depth -= 1
        current.append(char)
        if depth == 0 and current:
            block = "".join(current).strip()
            if block and block.startswith("["):
                matrices.append(block)
            current = []

    if len(matrices) < 4:
        return None, None, None, None, success

    def parse_matrix(text):
        """Parse [[a b c] [d e f] ...] into a list of lists."""
        # Remove outer brackets
        text = text.strip()
        if text.startswith("[["):
            text = text[1:-1]  # remove outer [ ]
        # Split into rows by ][
        rows = re.split(r'\]\s*\[', text)
        result = []
        for row in rows:
            row = row.strip().strip("[]")
            vals = [int(x) for x in row.split()]
            result.append(vals)
        return result

    U_data = parse_matrix(matrices[0])
    M_data = parse_matrix(matrices[1])
    T_data = parse_matrix(matrices[2])
    S_data = parse_matrix(matrices[3])

    U = Matrix(U_data)
    M = Matrix(M_data)
    T = Matrix(T_data)
    S_vec = [row[0] for row in S_data]

    return U, M, T, S_vec, success


def check_smith_form(S_vec):
    """Check that s_1 | s_2 | ... | s_n."""
    for i in range(len(S_vec) - 1):
        if S_vec[i + 1] % S_vec[i] != 0:
            return False
    return True


def check_integrality(A, U, M, T, S_vec):
    """Check that B = diag(A,I) * ... * diag(I, S^{-1}) is integral.

    B = [[A, A*M*S^{-1}], [U, (T + U*M)*S^{-1}]]

    So we need:
      (a) A * M * S^{-1} is integral  (column j divided by S[j])
      (b) (T + U*M) * S^{-1} is integral
    """
    n = A.rows
    AM = A * M
    UM = U * M
    TpUM = T + UM

    # Check A*M*S^{-1}: column j must be divisible by S[j]
    for j in range(n):
        sj = S_vec[j]
        if sj == 0:
            continue
        for i in range(n):
            if AM[i, j] % sj != 0:
                return False, f"A*M not divisible by S at ({i},{j}): {AM[i,j]} mod {sj}"

    # Check (T + U*M)*S^{-1}: column j must be divisible by S[j]
    for j in range(n):
        sj = S_vec[j]
        if sj == 0:
            continue
        for i in range(n):
            if TpUM[i, j] % sj != 0:
                return False, f"(T+U*M) not divisible by S at ({i},{j}): {TpUM[i,j]} mod {sj}"

    return True, ""


def check_unimodularity(A, U, M, T, S_vec):
    """Check |det B| = 1 by computing det B = det A / prod(S)."""
    det_A = det(A)
    prod_S = 1
    for s in S_vec:
        prod_S *= s
    if prod_S == 0:
        return False, "S has a zero diagonal entry"
    # det(B) * det(S) = det(A), so det(B) = det(A) / prod(S)
    if det_A % prod_S != 0:
        return False, f"det(A) = {det_A} not divisible by prod(S) = {prod_S}"
    det_B = det_A // prod_S
    if abs(det_B) != 1:
        return False, f"|det B| = {abs(det_B)}, expected 1"
    return True, ""


def run_test(n, binary):
    """Run all checks for a given prime n. Returns (pass, details)."""
    result = subprocess.run(
        [binary, str(n)],
        capture_output=True, text=True, timeout=300
    )
    if result.returncode != 0:
        return False, "smith-massager crashed"

    U, M, T, S_vec, success = parse_output(result.stdout, n)

    if not success:
        return False, "smith-massager reported failure"
    if U is None:
        return False, "could not parse output"

    details = []

    # Check 1: S is in Smith form
    if check_smith_form(S_vec):
        details.append("divisibility: ok")
    else:
        return False, "S is NOT in Smith form (divisibility fails)"

    # Check 2: S matches SymPy
    A = build_test_matrix(n)
    sympy_diag = [snf(A, domain=ZZ)[i, i] for i in range(n)]
    if list(sympy_diag) == list(S_vec):
        details.append("SymPy: ok")
    else:
        return False, f"S does not match SymPy: got {S_vec}, expected {list(sympy_diag)}"

    # Check 3: B is integral
    ok, msg = check_integrality(A, U, M, T, S_vec)
    if ok:
        details.append("integrality: ok")
    else:
        return False, f"B is NOT integral: {msg}"

    # Check 4: |det B| = 1
    ok, msg = check_unimodularity(A, U, M, T, S_vec)
    if ok:
        details.append("unimodularity: ok")
    else:
        return False, f"B is NOT unimodular: {msg}"

    return True, ", ".join(details)


def main():
    parser = argparse.ArgumentParser(description="Verify SmithMassagerC against SymPy")
    parser.add_argument("max_n", nargs="?", type=int, default=53,
                        help="Test all primes up to this value (default: 53)")
    parser.add_argument("--binary", default=DEFAULT_BINARY,
                        help="Path to smith-massager binary")
    args = parser.parse_args()

    primes = [p for p in PRIMES if p <= args.max_n]

    print(f"Verifying Smith massager for {len(primes)} primes up to {args.max_n}")
    print(f"Binary: {args.binary}")
    print(f"Checks: divisibility, SymPy, integrality of B, unimodularity of B")
    print()

    all_pass = True
    for n in primes:
        print(f"n={n:4d}: ", end="", flush=True)
        try:
            ok, details = run_test(n, args.binary)
        except subprocess.TimeoutExpired:
            ok, details = False, "TIMEOUT"
        except Exception as e:
            ok, details = False, f"ERROR: {e}"

        if ok:
            print(f"PASS  ({details})")
        else:
            print(f"FAIL  ({details})")
            all_pass = False

    print()
    if all_pass:
        print(f"All {len(primes)} tests passed.")
    else:
        print("Some tests FAILED.")

    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
