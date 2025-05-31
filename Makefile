.PHONY: default all clean lib

include defs.mk

default: lib examples
	
#all: lib shared examples tuning version
# Everything related to iherm under test/ are broken.
all: lib shared examples

##################

lib:
	mkdir -p lib
	mkdir -p bin
	mkdir -p objs
	$(MAKE) $(MAKEOPTS) --directory=src lib
	cp src/maple/extern.mpl lib/

shared:
	$(MAKE) $(MAKEOPTS) --directory=src shared

examples: lib shared
	$(MAKE) $(MAKEOPTS) --directory=examples

#tests: lib
#	$(MAKE) $(MAKEOPTS) --directory=tests tests

#version: lib
#	$(MAKE) $(MAKEOPTS) --directory=tests version

#tuning: lib
#	$(MAKE) $(MAKEOPTS) --directory=tests tuning

clean:
	$(MAKE) $(MAKEOPTS) --directory=src clean
	rm -f objs/* lib/* bin/*; true;
install:
	#mkdir $(INSTALL_DIR)
	cp -r lib/ $(INSTALL_DIR)

