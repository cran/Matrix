#-------------------------------------------------------------------------------
# UMFPACK makefile (for GNU make or original make)
#-------------------------------------------------------------------------------

# UMFPACK requires the AMD package to be in ../AMD

default: library

# compile all C code (except hb, fortran, and fortran64), including AMD and the
# MATLAB mexFunctions
all:
	( cd Source ; make )

# compile just the C-callable libraries and demo programs (not mexFunctions)
library:
	( cd Source ; make )

# remove object files, but keep the compiled programs and library archives
clean:
	( cd Source ; make clean )

# clean, and then remove compiled programs and library archives
purge:
	( cd Source ; make purge )

# create PDF documents for the original distribution
doc:
	#( cd ../AMD ; make doc )
	#( cd Doc    ; make )

distclean: purge

ccode: library

lib: library
