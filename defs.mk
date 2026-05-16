export BASEDIR := $(dir $(CURDIR)/$(word $(words $(MAKEFILE_LIST)),$(MAKEFILE_LIST)))

# Compiler.
export CC := cc

# Required libraries.  (Maple is optional; see MAPLEDIR below.)
export LDLIBS := -lopenblas -lgmp -lm -lflint

# Platform detection.
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Default install locations.  Override these on the command line if
# your libraries live elsewhere, e.g.
#   make GMP_LIB_DIR=/path/to/gmp/lib ...
ifeq ($(UNAME_S),Darwin)
  # macOS: assume Homebrew.
  ifeq ($(UNAME_M),arm64)
    HOMEBREW ?= /opt/homebrew/opt
  else
    HOMEBREW ?= /usr/local/opt
  endif
  GMP_LIB_DIR ?= $(HOMEBREW)/gmp/lib/
  GMP_INCLUDE_DIR ?= $(HOMEBREW)/gmp/include/
  OPENBLAS_LIB_DIR ?= $(HOMEBREW)/openblas/lib/
  OPENBLAS_INCLUDE_DIR ?= $(HOMEBREW)/openblas/include/
  FLINT_LIB_DIR ?= $(HOMEBREW)/flint/lib/
  FLINT_INCLUDE_DIR ?= $(HOMEBREW)/flint/include/flint/
  OS := APPLE
  ifeq ($(UNAME_M),arm64)
    ARCH := ARM64_MACOS
  else
    ARCH := UNIVERSAL_OSX
  endif
else
  # Linux: assume system packages (e.g. apt install libgmp-dev
  # libopenblas-dev libflint-dev).  Override LIB_DIR / INCLUDE_DIR
  # variables if your libraries live in a non-standard location.
  FLINT_INCLUDE_DIR ?= /usr/include/flint/
  OS := LINUX
  ifeq ($(UNAME_M),x86_64)
    ARCH := X86_64
  else ifeq ($(UNAME_M),arm64)
    ARCH := ARM64
  endif
endif

# Maple bindings (optional).  Set MAPLEDIR on the command line to
# enable, e.g.  make MAPLEDIR=/opt/maple2023/ shared
ifdef MAPLEDIR
  export MAPLEBIN := $(MAPLEDIR)/bin.$(if $(filter LINUX,$(OS)),$(ARCH)_$(OS),$(OS)_$(ARCH))/
  LDLIBS += -lmaplec
endif

export INSTALL_DIR := ./install/

export LDFLAGS := #empty
export CFLAGS  := #empty

export OBJDIR := $(BASEDIR)/objs/
export SRCDIR := $(BASEDIR)/src/

export STATICLIB := $(BASEDIR)/lib/libhnfproj.a
export SHAREDLIB := $(BASEDIR)/lib/libhnfproj.so
export SMITH := $(BASEDIR)/bin/smith-massager

ifdef OPENBLAS_LIB_DIR
  LDFLAGS += -L$(OPENBLAS_LIB_DIR)
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

CFLAGS += -fPIC
CFLAGS += -pedantic
CFLAGS += -Wno-all
CFLAGS += -Wno-extra
CFLAGS += -Wshadow
CFLAGS += -Wpointer-arith
CFLAGS += -Wcast-align
CFLAGS += -Wstrict-prototypes
CFLAGS += -Wmissing-prototypes
CFLAGS += -Wno-long-long
CFLAGS += -Wno-variadic-macros
CFLAGS += -Wno-implicit-function-declaration
CFLAGS += -Wno-int-conversion

ifdef NOTIMER
  CFLAGS += -DNOTIMER
endif
ifdef NOPRINT
  CFLAGS += -DNOPRINT
endif

ifdef DEBUG
  CFLAGS += -O0 -g -DDEBUG
else
  CFLAGS += -O3 -DNDEBUG
endif

ifdef THREAD
  CFLAGS += -DTHREAD
  LDLIBS += -lpthread
endif

ifeq ($(CC),icc)
  CFLAGS += -wd1782
endif

export MAKEOPTS := #empty
