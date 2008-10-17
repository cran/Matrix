#===============================================================================
# SuiteSparseQR/Lib/Makefile
#===============================================================================

include ${MkInclude}

default: lib

## C preprocessor: For C and C++ :
PKG_CPPFLAGS = -I../../UFconfig -I../../CHOLMOD/Include -I../Include -DNPARTITION

# no statement coverage, high optimization
# CFLAGS = -O3 \
# 	-Wall -W -Wshadow -Winline -Wno-unused-parameter \
# 	-Wredundant-decls -Wdisabled-optimization \
# 	-ansi -fexceptions

INC = ../Include/spqr.hpp ../Include/SuiteSparseQR_C.h \
	../Include/SuiteSparseQR_definitions.h \
	../Include/SuiteSparseQR.hpp Makefile

OBJ = \
    spqr_rmap.o \
    SuiteSparseQR_C.o \
    SuiteSparseQR_expert.o \
    spqr_parallel.o \
    spqr_kernel.o \
    spqr_analyze.o \
    spqr_assemble.o \
    spqr_cpack.o \
    spqr_csize.o \
    spqr_fcsize.o \
    spqr_debug.o \
    spqr_front.o \
    spqr_factorize.o \
    spqr_freenum.o \
    spqr_freesym.o \
    spqr_freefac.o \
    spqr_fsize.o \
    spqr_maxcolnorm.o \
    spqr_rconvert.o \
    spqr_rcount.o \
    spqr_rhpack.o \
    spqr_rsolve.o \
    spqr_stranspose1.o \
    spqr_stranspose2.o \
    spqr_hpinv.o \
    spqr_1fixed.o \
    spqr_1colamd.o \
    SuiteSparseQR.o \
    spqr_1factor.o \
    spqr_cumsum.o \
    spqr_shift.o \
    spqr_happly.o \
    spqr_panel.o \
    spqr_happly_work.o \
    SuiteSparseQR_qmult.o \
    spqr_trapezoidal.o \
    spqr_larftb.o \
    spqr_append.o \
    spqr_type.o \
    spqr_tol.o

LIB = ../../SPQR.a

lib: $(LIB)

$(LIB): $(OBJ)
	$(AR) -rucs $(LIB) $(OBJ)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f $(OBJ)

C = $(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS) $(SPQR_CONFIG)

spqr_1colamd.o: ../Source/spqr_1colamd.cpp
	$(C) -c ../Source/spqr_1colamd.cpp

spqr_1factor.o: ../Source/spqr_1factor.cpp
	$(C) -c ../Source/spqr_1factor.cpp

spqr_1fixed.o: ../Source/spqr_1fixed.cpp
	$(C) -c ../Source/spqr_1fixed.cpp

spqr_analyze.o: ../Source/spqr_analyze.cpp
	$(C) -c ../Source/spqr_analyze.cpp

spqr_parallel.o: ../Source/spqr_parallel.cpp
	$(C) -c ../Source/spqr_parallel.cpp

spqr_kernel.o: ../Source/spqr_kernel.cpp
	$(C) -c ../Source/spqr_kernel.cpp

spqr_append.o: ../Source/spqr_append.cpp
	$(C) -c ../Source/spqr_append.cpp

spqr_assemble.o: ../Source/spqr_assemble.cpp
	$(C) -c ../Source/spqr_assemble.cpp

spqr_cpack.o: ../Source/spqr_cpack.cpp
	$(C) -c ../Source/spqr_cpack.cpp

spqr_csize.o: ../Source/spqr_csize.cpp
	$(C) -c ../Source/spqr_csize.cpp

spqr_cumsum.o: ../Source/spqr_cumsum.cpp
	$(C) -c ../Source/spqr_cumsum.cpp

spqr_debug.o: ../Source/spqr_debug.cpp
	$(C) -c ../Source/spqr_debug.cpp

spqr_factorize.o: ../Source/spqr_factorize.cpp
	$(C) -c ../Source/spqr_factorize.cpp

spqr_fcsize.o: ../Source/spqr_fcsize.cpp
	$(C) -c ../Source/spqr_fcsize.cpp

spqr_freefac.o: ../Source/spqr_freefac.cpp
	$(C) -c ../Source/spqr_freefac.cpp

spqr_freenum.o: ../Source/spqr_freenum.cpp
	$(C) -c ../Source/spqr_freenum.cpp

spqr_freesym.o: ../Source/spqr_freesym.cpp
	$(C) -c ../Source/spqr_freesym.cpp

spqr_fsize.o: ../Source/spqr_fsize.cpp
	$(C) -c ../Source/spqr_fsize.cpp

spqr_happly.o: ../Source/spqr_happly.cpp
	$(C) -c ../Source/spqr_happly.cpp

spqr_panel.o: ../Source/spqr_panel.cpp
	$(C) -c ../Source/spqr_panel.cpp

spqr_happly_work.o: ../Source/spqr_happly_work.cpp
	$(C) -c ../Source/spqr_happly_work.cpp

spqr_hpinv.o: ../Source/spqr_hpinv.cpp
	$(C) -c ../Source/spqr_hpinv.cpp

spqr_larftb.o: ../Source/spqr_larftb.cpp
	$(C) -c ../Source/spqr_larftb.cpp

spqr_rconvert.o: ../Source/spqr_rconvert.cpp
	$(C) -c ../Source/spqr_rconvert.cpp

spqr_rcount.o: ../Source/spqr_rcount.cpp
	$(C) -c ../Source/spqr_rcount.cpp

spqr_rhpack.o: ../Source/spqr_rhpack.cpp
	$(C) -c ../Source/spqr_rhpack.cpp

spqr_rsolve.o: ../Source/spqr_rsolve.cpp
	$(C) -c ../Source/spqr_rsolve.cpp

spqr_shift.o: ../Source/spqr_shift.cpp
	$(C) -c ../Source/spqr_shift.cpp

spqr_stranspose1.o: ../Source/spqr_stranspose1.cpp
	$(C) -c ../Source/spqr_stranspose1.cpp

spqr_stranspose2.o: ../Source/spqr_stranspose2.cpp
	$(C) -c ../Source/spqr_stranspose2.cpp

spqr_trapezoidal.o: ../Source/spqr_trapezoidal.cpp
	$(C) -c ../Source/spqr_trapezoidal.cpp

spqr_type.o: ../Source/spqr_type.cpp
	$(C) -c ../Source/spqr_type.cpp

spqr_front.o: ../Source/spqr_front.cpp
	$(C) -c ../Source/spqr_front.cpp

SuiteSparseQR_expert.o: ../Source/SuiteSparseQR_expert.cpp
	$(C) -c ../Source/SuiteSparseQR_expert.cpp

spqr_maxcolnorm.o: ../Source/spqr_maxcolnorm.cpp
	$(C) -c ../Source/spqr_maxcolnorm.cpp

SuiteSparseQR_qmult.o: ../Source/SuiteSparseQR_qmult.cpp
	$(C) -c ../Source/SuiteSparseQR_qmult.cpp

SuiteSparseQR.o: ../Source/SuiteSparseQR.cpp
	$(C) -c ../Source/SuiteSparseQR.cpp

spqr_tol.o: ../Source/spqr_tol.cpp
	$(C) -c ../Source/spqr_tol.cpp

SuiteSparseQR_C.o: ../Source/SuiteSparseQR_C.cpp
	$(C) -c ../Source/SuiteSparseQR_C.cpp

spqr_rmap.o: ../Source/spqr_rmap.cpp
	$(C) -c ../Source/spqr_rmap.cpp
