# default: ccolamd_example ccolamd_l_example libccolamd.a
default: libccolamd.a

include ../UFconfig/UFconfig.mk

purge: distclean

distclean: clean2
	- $(RM) libccolamd.a

clean2: clean

# creates libccolamd.a, a C-callable CCOLAMD library
libccolamd.a:  ccolamd.c ccolamd_global.c ccolamd.h
	$(CC) $(CFLAGS) ${CPICFLAGS} -c ccolamd_global.c
	$(CC) $(CFLAGS) ${CPICFLAGS} -c ccolamd.c
	$(AR) libccolamd.a ccolamd.o ccolamd_global.o
	- cp libccolamd.a ../CCOLAMD.a

ccode: libccolamd.a

library: libccolamd.a

lib: libccolamd.a
