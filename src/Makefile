#
pkg = Matrix

PKG_CFLAGS = -I./Metis
PKG_LIBS = ${LAPACK_LIBS} ${BLAS_LIBS}
include SOURCES_C.mkf
OBJECTS = $(SOURCES_C:.c=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = Metis
SUBLIBS = $(SUBDIRS:=.a)
SUBSTAMP = $(SUBDIRS:=.stamp)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): $(SUBSTAMP)
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) lib) || exit 1; \
	  touch $${d}.stamp; \
	done

$(SUBSTAMP): 
	touch $(SUBSTAMP)

clean:
	@-rm -rf .libs _libs
	@-rm -f *.o *.so *.a
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) clean) || exit 1; \
	  rm -f $${d}.stamp; \
	done

include DEPS.mkf
