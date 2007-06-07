#===============================================================================
# CHOLOMD/Lib/Makefile: for compiling the CHOLMOD library
#===============================================================================
include ${MkInclude}

default: lib

PKG_CFLAGS = -I../../AMD/Include -I../../AMD/Source \
        -I../../COLAMD/Include \
	-I../Include -I../../UFconfig -DNPARTITION
#	 -I../../CCOLAMD \
#        -I../../CAMD/Include -I../../CAMD/Source \
#	-I../../Metis \

#-------------------------------------------------------------------------------
# ../Include/ directory contains all include files:
#-------------------------------------------------------------------------------

INC =   ../Include/cholmod.h \
	../Include/cholmod_blas.h \
	../Include/cholmod_check.h \
	../Include/cholmod_cholesky.h \
	../Include/cholmod_complexity.h \
	../Include/cholmod_config.h \
	../Include/cholmod_core.h \
	../Include/cholmod_internal.h \
	../Include/cholmod_matrixops.h \
	../Include/cholmod_modify.h \
	../Include/cholmod_supernodal.h \
	../Include/cholmod_template.h
#	../Include/cholmod_partition.h \

#-------------------------------------------------------------------------------
# The 7 CHOLMOD library modules (int, double)
#-------------------------------------------------------------------------------

CORE = cholmod_aat.o cholmod_add.o cholmod_band.o \
	cholmod_change_factor.o cholmod_common.o cholmod_complex.o \
	cholmod_copy.o cholmod_dense.o cholmod_error.o cholmod_factor.o \
	cholmod_memory.o cholmod_sparse.o \
	cholmod_transpose.o cholmod_triplet.o

CHECK = cholmod_check.o cholmod_read.o cholmod_write.o

CHOLESKY = cholmod_amd.o cholmod_analyze.o cholmod_colamd.o \
	cholmod_etree.o cholmod_factorize.o cholmod_postorder.o \
	cholmod_rcond.o cholmod_resymbol.o cholmod_rowcolcounts.o \
	cholmod_rowfac.o cholmod_solve.o cholmod_spsolve.o

MATRIXOPS = cholmod_drop.o cholmod_horzcat.o cholmod_norm.o \
	cholmod_scale.o cholmod_sdmult.o cholmod_ssmult.o \
	cholmod_submatrix.o cholmod_vertcat.o cholmod_symmetry.o

PARTITION = cholmod_ccolamd.o cholmod_csymamd.o \
	cholmod_metis.o cholmod_nesdis.o cholmod_camd.o

MODIFY = cholmod_rowadd.o cholmod_rowdel.o cholmod_updown.o

SUPERNODAL = cholmod_super_numeric.o cholmod_super_solve.o \
	cholmod_super_symbolic.o

OBJS = $(CORE) $(CHOLESKY) $(MATRIXOPS) $(MODIFY) $(SUPERNODAL) $(CHECK) #$(PARTITION) 

LIB = ../../CHOLMOD.a

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rucs $(LIB) $(OBJS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f $(OBJS) 

$(OBJ): $(INC)

#-------------------------------------------------------------------------------
# Check Module:
#-------------------------------------------------------------------------------

cholmod_check.o: ../Check/cholmod_check.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_read.o: ../Check/cholmod_read.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_write.o: ../Check/cholmod_write.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# Core Module:
#-------------------------------------------------------------------------------

cholmod_common.o: ../Core/cholmod_common.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_dense.o: ../Core/cholmod_dense.c ../Core/t_cholmod_dense.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_factor.o: ../Core/cholmod_factor.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_change_factor.o: ../Core/cholmod_change_factor.c \
	../Core/t_cholmod_change_factor.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_memory.o: ../Core/cholmod_memory.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_sparse.o: ../Core/cholmod_sparse.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_complex.o: ../Core/cholmod_complex.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_transpose.o: ../Core/cholmod_transpose.c ../Core/t_cholmod_transpose.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_band.o: ../Core/cholmod_band.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_copy.o: ../Core/cholmod_copy.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_triplet.o: ../Core/cholmod_triplet.c ../Core/t_cholmod_triplet.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_error.o: ../Core/cholmod_error.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_aat.o: ../Core/cholmod_aat.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_add.o: ../Core/cholmod_add.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# Cholesky Module:
#-------------------------------------------------------------------------------

cholmod_amd.o: ../Cholesky/cholmod_amd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_analyze.o: ../Cholesky/cholmod_analyze.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_colamd.o: ../Cholesky/cholmod_colamd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_etree.o: ../Cholesky/cholmod_etree.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_factorize.o: ../Cholesky/cholmod_factorize.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_postorder.o: ../Cholesky/cholmod_postorder.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_rcond.o: ../Cholesky/cholmod_rcond.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_resymbol.o: ../Cholesky/cholmod_resymbol.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_rowcolcounts.o: ../Cholesky/cholmod_rowcolcounts.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_solve.o: ../Cholesky/cholmod_solve.c ../Cholesky/t_cholmod_lsolve.c \
	../Cholesky/t_cholmod_ltsolve.c ../Cholesky/t_cholmod_solve.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_spsolve.o: ../Cholesky/cholmod_spsolve.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_rowfac.o: ../Cholesky/cholmod_rowfac.c ../Cholesky/t_cholmod_rowfac.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# Partition Module:
#-------------------------------------------------------------------------------

cholmod_ccolamd.o: ../Partition/cholmod_ccolamd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_csymamd.o: ../Partition/cholmod_csymamd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_metis.o: ../Partition/cholmod_metis.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_nesdis.o: ../Partition/cholmod_nesdis.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_camd.o: ../Partition/cholmod_camd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# MatrixOps Module:
#-------------------------------------------------------------------------------

cholmod_horzcat.o: ../MatrixOps/cholmod_horzcat.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_norm.o: ../MatrixOps/cholmod_norm.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_scale.o: ../MatrixOps/cholmod_scale.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_drop.o: ../MatrixOps/cholmod_drop.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_sdmult.o: ../MatrixOps/cholmod_sdmult.c \
	../MatrixOps/t_cholmod_sdmult.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_ssmult.o: ../MatrixOps/cholmod_ssmult.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_submatrix.o: ../MatrixOps/cholmod_submatrix.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_vertcat.o: ../MatrixOps/cholmod_vertcat.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_symmetry.o: ../MatrixOps/cholmod_symmetry.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# Modify Module:
#-------------------------------------------------------------------------------

cholmod_rowadd.o: ../Modify/cholmod_rowadd.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_rowdel.o: ../Modify/cholmod_rowdel.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_updown.o: ../Modify/cholmod_updown.c \
	../Modify/t_cholmod_updown.c ../Modify/t_cholmod_updown_numkr.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# Supernodal Module:
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

cholmod_super_numeric.o: ../Supernodal/cholmod_super_numeric.c \
	../Supernodal/t_cholmod_super_numeric.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_super_symbolic.o: ../Supernodal/cholmod_super_symbolic.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
cholmod_super_solve.o: ../Supernodal/cholmod_super_solve.c \
	../Supernodal/t_cholmod_super_solve.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@

