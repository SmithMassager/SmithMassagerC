# SmithMassagerC

A C implementation of the Las Vegas algorithm for computing the
Smith normal form of a nonsingular integer matrix, based on:

> S. Birmpilis, G. Labahn, A. Storjohann.
> *A Las Vegas algorithm for computing the Smith form of a nonsingular
>  integer matrix.*  ISSAC 2020.

The algorithm runs in time softly equivalent to integer matrix
multiplication: `O(n^ω · B(log n + log ‖A‖) · (log n)²)` bit operations.

## Dependencies

Required:

- [GMP](https://gmplib.org/) — arbitrary precision arithmetic
- [FLINT](https://flintlib.org/) — integer matrix operations (v3 or later)
- [OpenBLAS](https://www.openblas.net/) — BLAS for double-precision matrix multiply

Optional:

- [Maple](https://www.maplesoft.com/) — only needed for the Maple bindings
  under `src/maple/`

The Integer Matrix Library (IML) is vendored under `src/iml-1.0.3/`.

## Quick start

### macOS (Apple Silicon, Homebrew defaults)

```
brew install gmp openblas flint
make lib
./bin/smith-massager 53
```

### macOS (Intel)

Edit `defs.mk` and set `HOMEBREW := /usr/local/opt`, then:

```
make lib
./bin/smith-massager 53
```

### Linux (Ubuntu / Debian)

```
sudo apt install libgmp-dev libopenblas-dev libflint-dev
make lib
./bin/smith-massager 53
```

The `smith-massager` binary takes a prime `n` as argument and computes
the Smith massager of the n × n test matrix `A[i,j] = i^j mod n`.
Output includes the matrices `(U, M, T, S)` and a `success: 1` line.

## Build targets

| Target | Description |
|---|---|
| `make lib` | Build static library `lib/libhnfproj.a` and `bin/smith-massager` |
| `make shared` | Build shared library `lib/libhnfproj.so` |
| `make examples` | Build example binaries under `examples/` |
| `make clean` | Remove `objs/`, `lib/`, `bin/` |
| `make install` | Install to `INSTALL_DIR` (see `defs.mk`) |

## Configuration

All build configuration is in `defs.mk`.  Key settings:

- `GMP_LIB_DIR`, `GMP_INCLUDE_DIR` — GMP install location
- `OPENBLAS_LIB_DIR`, `OPENBLAS_INCLUDE_DIR` — OpenBLAS install location
- `FLINT_LIB_DIR`, `FLINT_INCLUDE_DIR` — FLINT install location
- `MAPLEDIR` — Maple install location (set to enable Maple bindings)
- `THREAD=1` — enable pthread parallelism for RNS operations
- `DEBUG=1` — build with `-O0 -g` and debug prints

## Source layout

```
src/core/         Smith massager algorithm + RNS infrastructure
src/lift/         Double-plus-one p-adic lifting
src/spinv/        Sparse inverse, high-order residues, unimodularity certification
src/iherm/        Hermite normal form (independent module)
src/iml-1.0.3/    Vendored IML library
src/maple/        Optional Maple bindings
examples/         Example programs
tests/            Test scripts
```

## License

BSD 3-clause.  See [LICENSE](LICENSE).
