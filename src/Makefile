#
pkg = Matrix

PKG_CFLAGS = -I./UFconfig
## we use the BLAS and the LAPACK library:
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

ifdef build_under_Windows
  DLLLIBS = $(PKG_LIBS) -lR
  SHLIB_EXT = .dll
  PKG_HOME = $(shell $(RHOME)/src/gnuwin32/Rpwd.exe .)
  MkInclude = $(PKG_HOME)/Win.mk
else
  MkInclude = $(R_HOME)/etc${R_ARCH}/Makeconf
endif

include SOURCES_C.mkf

OBJECTS = $(SOURCES_C:.c=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD COLAMD AMD
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
ifdef build_under_Windows
	$(DLL) --shared -s -o $@ $(OBJECTS) $(SUBLIBS) $(DLLLIBS)
else
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)
endif
all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude=$(MkInclude) $(MAKE) lib) || exit 1; \
	  touch $${d}.stamp; \
	done

$(SUBSTAMP):
	touch $(SUBSTAMP)

clean:
	@-rm -rf .libs _libs
	@-rm -f *.o $(SHLIB) *.a $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) clean) || exit 1; \
	  rm -f $${d}.stamp; \
	done

include DEPS.mkf
