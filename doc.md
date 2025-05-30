<h1 align="center">  Documentation  </h1>

> Requires gcc compiler, openblas library, GMP, maple (optional but needed for
> the Smith Massager project).

### Installation & Usage

```
> make CC=gcc OPENBLAS_LIB_DIR=your_blas_lib OPENBLAS_INCLUDE_DIR=your_blas_header GMP_LIB_DIR=your_gmp_lib GMP_INCLUDE_DIR=your_gmp_header MAPLEDIR=your_maple_dir
> make install
```

Usually Openblas and GMP are already installed on your computer. Some common
places to find it:
```
> cd /lib/
> find . -name '*gmp*'
...
./x86_64-linux-gnu/libgmp.so.10
...
### Similarly for blas. And so your GMP_LIB_DIR=/lib/x86_64-linux-gnu/libgmp.so.10
> cd /usr
> find . -name '*gmp*'
...
./local/include/gmp.h
...
### Similarly for blas. And so your GMP_INCLUDE_DIR=/usr/local/include/gmp.h
```

The last ``MAPLEDIR`` is optional to enable maple interfaces. And this can
directory can be directly found by entering the following command ``which
maple`` a sample output:

```
> which maple
/opt/maple2023/bin/maple
```

So ``MAPLEDIR=/opt/maple2023/``.

**NOTE** If you have get a compilation ever complaining about the file
``maple_c.h`` or ``maplec.h`` is not found. Then please change the ``#include
"maplec.h"`` to ``#include "maple_c.h"`` inside ``include/maple_call.h`` or
vice versa. This issue is because at some point Maple decided to change their
filename.

**NOTE** the default install directory should be ``../maple/`` which will
automatically install the library to the Smith Massager project and that
project will need maple support.

Main features:
- ``mpzMatrix_t``, matrix of gmp integers allowing large entries. And basic matrix
  operations.
- ``rnsMatrix_t``, residue number system matrix, for each residue number the corresponding matrix is
  using built in data types that takes advantage of BLAS library. There are
  interfaces that converts between ``mpzMatrix_t`` and ``rnsMatrix_t``. There
  is a Chinese Remainder function to compute the inverse of a integer matrix.
- DoublePlusOneLift, implementation of the double plus one lifting which are
  used to determine Unimodular matrices and High order liftings.
- Maple Interfaces, Maple interface that uses this project for fast integer
  certificate, unimodular certificate, iml solve, iml matrix multiplication.
