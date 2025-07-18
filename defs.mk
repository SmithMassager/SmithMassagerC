export BASEDIR := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

export CC := gcc
export LDLIBS  := -lopenblas -lgmp -lm -lmaplec -lflint
export prefix := /Users/ziwenwang/repo/flint-nix-builds/
export OPENBLAS := $(prefix)/openblas/
export OPENBLAS_LIB_DIR := $(OPENBLAS)/lib/
export OPENBLAS_INCLUDE_DIR := $(OPENBLAS)/include/
export GMP := $(prefix)/gmp/
export GMP_LIB_DIR := $(GMP)/lib/
export GMP_INCLUDE_DIR := $(GMP)/include/
export FLINT := $(prefix)/flint/
export FLINT_LIB_DIR := $(FLINT)/lib/
export FLINT_INCLUDE_DIR := $(FLINT)/include/

export INSTALL_DIR := ../maple/

#export THREAD := true

export MAPLEDIR := /Library/Frameworks/Maple.framework/Versions/2025/

export LDFLAGS := #empty
export CFLAGS  := #empty

export OBJDIR := $(BASEDIR)objs/
export SRCDIR := $(BASEDIR)src/

export LIFTLIB := $(BASEDIR)lib/libhnfproj.a
export SHAREDLIB := $(BASEDIR)lib/libhnfproj.so

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Linux)
  OS := LINUX
  ifeq ($(UNAME_M),x86_64)
    ARCH := X86_64
  else ifeq ($(UNAME_M),arm64)
    ARCH := ARM64
	endif
  MAPLEBIN := $(MAPLEDIR)/bin.$(ARCH)_$(OS)/
else ifeq ($(UNAME_S),Darwin)
  OS := APPLE
  ifeq ($(UNAME_M),x86_64)
    ARCH := UNIVERSAL_OSX
  else ifeq ($(UNAME_M),arm64)
    ARCH := ARM64_MACOS
	endif
  MAPLEBIN := $(MAPLEDIR)/bin.$(OS)_$(ARCH)/
else
    OS := unknown
    ARCH := unkonwn
endif

CFLAGS += -I$(MAPLEDIR)/extern/include/ -Wl,-rpath,$(MAPLEBIN) -Wl,-undefined,dynamic_lookup

LDFLAGS += -L$(MAPLEDIR)/lib/ -L$(MAPLEBIN)

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
CFLAGS  += -Wno-all
CFLAGS  += -Wno-extra
CFLAGS  += -Wshadow
CFLAGS  += -Wpointer-arith
CFLAGS  += -Wcast-align
CFLAGS  += -Wstrict-prototypes
CFLAGS  += -Wmissing-prototypes
CFLAGS  += -Wno-long-long
CFLAGS  += -Wno-variadic-macros
CFLAGS  += -Wno-implicit-function-declaration
CFLAGS  += -Wno-int-conversion

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
