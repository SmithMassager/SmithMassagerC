.PHONY: core lift maple spinv iherm clean

CFLAGS += -I$(SRCDIR)iml-1.0.3/src/ -I$(SRCDIR)/include/ -I$(SRCDIR)/lift/ -I$(SRCDIR)/core/ -I$(SRCDIR)/spinv/

lib: $(SMITH) $(SHAREDLIB) $(STATIC)

shared: $(SHAREDLIB)
static: $(STATIC)
smith: $(SMITH)

BASE_DEPS = iml core lift spinv iherm
MAPLE_DEPS =

$(OBJDIR)*.o: $(BASE_DEPS)

$(SMITH): $(BASE_DEPS)
		$(CC) $(OBJDIR)/*.o $(CFLAGS) $(LDFLAGS) $(LDLIBS) -o $@

$(STATIC): $(BASE_DEPS)
		$(AR) -r -c -s $@ $(OBJDIR)/*.o

ifdef MAPLEDIR
$(SHAREDLIB): CFLAGS += -I$(MAPLEDIR)/extern/include/ -Wl,-rpath,$(MAPLEBIN) -Wl,-undefined,dynamic_lookup
$(SHAREDLIB): LDFLAGS += -L$(MAPLEDIR)/lib/ -L$(MAPLEBIN)
$(SHAREDLIB): LDLIBS += -lmaplec  # Or -lmaple if that's your lib; confirm from Maple docs
MAPLE_DEPS = maple
endif

$(SHAREDLIB): $(BASE_DEPS) $(MAPLE_DEPS)
		$(CC) $(OBJDIR)/*.o $(CFLAGS) $(MAPLEOPTS) -shared $(LDFLAGS) $(LDLIBS) -o $@

iml:
	$(MAKE) $(MAKEOPTS) --directory=iml-1.0.3

core:
	$(MAKE) $(MAKEOPTS) --directory=core

lift:
	$(MAKE) $(MAKEOPTS) --directory=lift

spinv:
	$(MAKE) $(MAKEOPTS) --directory=spinv

iherm:
	$(MAKE) $(MAKEOPTS) --directory=iherm

maple:
	$(MAKE) $(MAKEOPTS) --directory=maple

clean:
	$(MAKE) $(MAKEOPTS) --directory=iml-1.0.3 clean
