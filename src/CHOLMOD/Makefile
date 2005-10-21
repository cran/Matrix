#-------------------------------------------------------------------------------
# CHOLMOD Makefile
#-------------------------------------------------------------------------------

default: all

# Compile the C-callable libraries and the Demo programs.
all:
	( cd Lib ; make )

# Compile the C-callable libraries only.
library:
	( cd Lib ; make )

lib: library

# Remove all files not in the original distribution
purge:
	( cd Lib ; make purge )

# Remove all files not in the original distribution, except keep the 
# compiled libraries.
clean:
	( cd Lib ; make clean )

distclean: purge

ccode: all

