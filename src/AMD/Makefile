#------------------------------------------------------------------------------
# AMD Makefile (for GNU Make or original make)
#------------------------------------------------------------------------------

default: library

# Compile all C code, including the C-callable routine and the mexFunctions.
# Do not compile the FORTRAN versions.
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

distclean: purge

# create PDF documents for the original distribution
doc:
	#( cd Doc    ; make )

# get ready for distribution
dist: purge
	#( cd Doc    ; make )

ccode: library

lib: library

# compile the MATLAB mexFunction
matlab:
	( cd MATLAB ; make )
