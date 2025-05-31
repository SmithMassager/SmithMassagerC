export BASEDIR := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

export CC      := gcc
export LDLIBS  := -lopenblas -lgmp -lm -lmaplec -lflint
export OPENBLAS := $(HOME)/OpenBLAS/install/
export OPENBLAS_LIB_DIR := $(OPENBLAS)/lib/
export OPENBLAS_INCLUDE_DIR := $(OPENBLAS)/include/
export GMP := $(HOME)/gmp-6.3.0/install/
export GMP_LIB_DIR := $(GMP)/lib/
export GMP_INCLUDE_DIR := $(GMP)/include/
export FLINT := $(HOME)/flint-3.2.2/install/
export FLINT_LIB_DIR := $(FLINT)/lib/
export FLINT_INCLUDE_DIR := $(FLINT)/include/flint/

export INSTALL_DIR := ../maple/

#export THREAD := true

export MAPLEDIR := /opt/maple2023/

export LDFLAGS := #empty
export CFLAGS  := #empty

export OBJDIR := $(BASEDIR)objs/
export SRCDIR := $(BASEDIR)src/

export LIFTLIB := $(BASEDIR)lib/libhnfproj.a
export SHAREDLIB := $(BASEDIR)lib/libhnfproj.so


CFLAGS += -I$(MAPLEDIR)/extern/include/ -Wl,-rpath=$(MAPLEDIR)/bin.X86_64_LINUX/

LDFLAGS += -L$(MAPLEDIR)/lib/ -L$(MAPLEDIR)/bin.X86_64_LINUX/

ifdef OPENBLAS_LIB_DIR
  LDFLAGS += -L$(OPENBLAS_LIB_DIR) -L$(HOME)/lib/
endif
ifdef OPENBLAS_INCLUDE_DIR
  CFLAGS += -I$(OPENBLAS_INCLUDE_DIR)
endif
ifdef GMP_LIB_DIR
  LDFLAGS += -L$(GMP_LIB_DIR)
endif
ifdef GMP_INCLUDE_DIR
  CFLAGS += -I$(GMP_INCLUDE_DIR)
endif
ifdef FLINT_LIB_DIR
  LDFLAGS += -L$(FLINT_LIB_DIR)
endif
ifdef FLINT_INCLUDE_DIR
  CFLAGS += -I$(FLINT_INCLUDE_DIR)
endif


CFLAGS  += -fPIC
CFLAGS  += -pedantic
CFLAGS  += -Wall
CFLAGS  += -Wextra
CFLAGS  += -Wshadow
CFLAGS  += -Wpointer-arith
CFLAGS  += -Wcast-align
CFLAGS  += -Wstrict-prototypes
CFLAGS  += -Wmissing-prototypes
CFLAGS  += -Wno-long-long
CFLAGS  += -Wno-variadic-macros

ifdef NOTIMER
  CFLAGS  += -DNOTIMER
endif
ifdef NOPRINT
  CFLAGS  += -DNOPRINT
endif

ifdef DEBUG
  CFLAGS  += -O0
  CFLAGS  += -g
  CFLAGS  += -DDEBUG
else
  CFLAGS  += -O3
  CFLAGS  += -DNDEBUG
endif

ifdef THREAD
  CFLAGS += -DTHREAD
  LDLIBS += -lpthread
endif

ifeq ($(CC),icc)
  CFLAGS += -wd1782 #pragma once is okay
endif

#CFLAGS += -DHAVE_VERBOSE_MODE -DHAVE_TIME_H

export MAKEOPTS := #empty
