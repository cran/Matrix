#
pkg = Matrix

PKG_CFLAGS = -I./Metis
PKG_LIBS = ${LAPACK_LIBS} ${BLAS_LIBS}
SOURCES_C =  \
	    LU.c \
            Mutils.c \
	    Metis_utils.c \
	    cscBlocked.c \
            cscMatrix.c \
            dense.c \
            factorizations.c \
            geMatrix.c \
	    init.c \
	    ldl.c \
	    lmeRep.c \
            poMatrix.c \
            sscCrosstab.c \
            sscMatrix.c \
            ssclme.c \
            syMatrix.c \
            trMatrix.c \
            triplet.c \
            triplet_to_col.c \
            tscMatrix.c \
            # flame.c
#SOURCES_F = foo.f
#OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o)
OBJECTS = $(SOURCES_C:.c=.o)
SHLIB = $(pkg)$(SHLIB_EXT)
SUBDIRS = Metis
SUBLIBS = $(SUBDIRS:=.a)

$(SHLIB): $(OBJECTS) $(SUBLIBS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(SUBLIBS) $(ALL_LIBS)

all: $(SHLIB)

mostlyclean: clean

$(SUBLIBS): 
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) lib) || exit 1; \
	done

clean:
	@-rm -rf .libs _libs
	@-rm -f *.o

LU.o: LU.c LU.h trMatrix.h Mutils.h
Metis_utils.o: Metis_utils.c Metis_utils.h \
  Metis/metis.h Metis/defs.h \
  Metis/struct.h Metis/macros.h Metis/rename.h Metis/proto.h
Mutils.o: Mutils.c Mutils.h triplet_to_col.h
cscMatrix.o: cscMatrix.c cscMatrix.h Mutils.h
dense.o: dense.c dense.h
factorizations.o: factorizations.c factorizations.h Mutils.h
flame.o: flame.c flame.h
geMatrix.o: geMatrix.c geMatrix.h 
geMutils.o: geMutils.c geMutils.h
ldl.o: ldl.c ldl.h
pdDiag.o: pdDiag.c Mutils.h
pdIdent.o: pdIdent.c Mutils.h
pdLogChol.o: pdLogChol.c Mutils.h
pdMat.o: pdMat.c Mutils.h
pdNatural.o: pdNatural.c Mutils.h
poMatrix.o: poMatrix.c poMatrix.h
sscChol.o: sscChol.c sscChol.h tscMatrix.h Mutils.h
sscCrosstab.o: sscCrosstab.c sscCrosstab.h Mutils.h
sscMatrix.o: sscMatrix.c sscMatrix.h 
ssclme.o: ssclme.c ssclme.h sscCrosstab.h Mutils.h ldl.h
syMatrix.o: syMatrix.c syMatrix.h geMatrix.h 
trMatrix.o: trMatrix.c trMatrix.h
triplet.o: triplet.c triplet.h Mutils.h
tscMatrix.o: tscMatrix.c tscMatrix.h Mutils.h
