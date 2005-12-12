#
pkg = Matrix

PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS)
include SOURCES_C.mkf

## zherk.f and lsame.f are part of BLAS, but not of R's BLAS sources
ifeq (, $(strip $(BLAS_LIBS)))
SOURCES_BLAS = zherk.f lsame.f
else
SOURCES_BLAS =
endif

## zpotrf.f and zpotf2.f are part of LAPACK, but not of R's LAPACK sources
ifeq (, $(findstring -lRlapack, $(LAPACK_LIBS)))
SOURCES_LAPACK =
else
SOURCES_LAPACK = zpotf2.f zpotrf.f 
endif

SOURCES_F = $(SOURCES_BLAS) $(SOURCES_LAPACK)

OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD UMFPACK COLAMD CCOLAMD AMD Metis LDL
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude=$(R_HOME)/etc/Makeconf $(MAKE) lib) || exit 1; \
	  touch $${d}.stamp; \
	done

$(SUBSTAMP):
	touch $(SUBSTAMP)

clean:
	@-rm -rf .libs _libs
	@-rm -f *.o *.so *.a $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) clean) || exit 1; \
	  rm -f $${d}.stamp; \
	done

include DEPS.mkf
