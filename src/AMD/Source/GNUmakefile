#-------------------------------------------------------------------------------
# AMD Makefile for compiling on Unix systems (for GNU make only)
#-------------------------------------------------------------------------------

default: ../Lib/libamd.a

include ../../UFconfig/UFconfig.mk

C = $(CC) $(CFLAGS) ${CPICFLAGS} $(CONFIG) -I../Include

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = ../Include/amd.h amd_internal.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
#AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_global.o: amd_global.c $(INC)
	$(C) -c $< -o $@

amd_i_%.o: amd_%.c $(INC)
	$(C) -DDINT -c $< -o $@

amd_l_%.o: amd_%.c $(INC)
	$(C) -DDLONG -c $< -o $@

#-------------------------------------------------------------------------------
# Create the libamd.a library (C versions only)
#-------------------------------------------------------------------------------

../Lib/libamd.a: amd_global.o $(AMDI) # $(AMDL)
	$(AR) ../Lib/libamd.a $^
	- $(RANLIB) ../Lib/libamd.a
	- cp ../Lib/libamd.a ../../AMD.a

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

purge: clean
	- $(RM) ../Lib/libamd.a ../Lib/libamdf77.a
