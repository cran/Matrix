#
pkg = Matrix

## we use the BLAS and the LAPACK library:
PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)
include SOURCES_C.mkf

## zpotrf.f and zpotf2.f are part of LAPACK, but not of R's LAPACK sources
## get rid of this, once we have 'Depends: R (>= 2.5.0)':
##ifeq (, $(findstring -lRlapack, $(LAPACK_LIBS)))
##SOURCES_LAPACK =
##else
##SOURCES_LAPACK = zpotf2.f zpotrf.f zlacgv.f
##endif

##SOURCES_F = $(SOURCES_LAPACK)

##OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
OBJECTS = $(SOURCES_C:.c=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = CHOLMOD COLAMD AMD
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && MkInclude=$(R_HOME)/etc${R_ARCH}/Makeconf $(MAKE) lib) || exit 1; \
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
