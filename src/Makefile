#
pkg = Matrix

PKG_CFLAGS = -DUSE_CHOLMOD
PKG_LIBS = ${LAPACK_LIBS} ${BLAS_LIBS}
include SOURCES_C.mkf

## zherk.f is part of BLAS, but not of R's BLAS sources
SOURCES_F = zpotf2.f zpotrf.f zherk.f

OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD UMFPACK Metis AMD COLAMD LDL CCOLAMD
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

UFconfig.stamp:
	cp UFconfig/UFconfig_unix.mk UFconfig/UFconfig.mk
	touch UFconfig.stamp

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): UFconfig.stamp $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) lib) || exit 1; \
	  touch $${d}.stamp; \
	done

$(SUBSTAMP):
	touch $(SUBSTAMP)

clean: UFconfig.stamp
	@-rm -rf .libs _libs
	@-rm -f *.o *.so *.a $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) clean) || exit 1; \
	  rm -f $${d}.stamp; \
	done

include DEPS.mkf
