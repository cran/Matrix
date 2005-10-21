#default: libcolamd.a colamd_example colamd_l_example
default: libcolamd.a

include ../UFconfig/UFconfig.mk

colamd_example: colamd_example.c libcolamd.a
	$(CC) $(CFLAGS) -o colamd_example colamd_example.c libcolamd.a -lm
	- ./colamd_example > my_colamd_example.out
	- diff colamd_example.out my_colamd_example.out

colamd_l_example: colamd_l_example.c libcolamd.a
	$(CC) $(CFLAGS) -o colamd_l_example colamd_l_example.c libcolamd.a -lm
	- ./colamd_l_example > my_colamd_l_example.out
	- diff colamd_example.out my_colamd_example.out

purge: distclean

distclean: clean2
	- $(RM) libcolamd.a

clean2: clean
	- $(RM) *.o *.dll colamd_example colamd_l_example
	- $(RM) colamdmex.mex* symamdmex.mex*
	- $(RM) colamdtestmex.mex* symamdtestmex.mex*
	- $(RM) my_colamd_example.out my_colamd_l_example.out

# creates libcolamd.a, a C-callable COLAMD library
libcolamd.a:  colamd.c colamd_global.c colamd.h
	$(CC) $(CFLAGS) ${CPICFLAGS} -c colamd_global.c
	$(CC) $(CFLAGS) ${CPICFLAGS} -c colamd.c
	$(AR) libcolamd.a colamd.o colamd_global.o
	- cp libcolamd.a ../COLAMD.a

ccode: libcolamd.a

library: libcolamd.a

lib: library
