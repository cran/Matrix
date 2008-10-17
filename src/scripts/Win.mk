# -*- Makefile -*-
## Note: This is included as 'MkInclude' from several */Makefile s

include $(RHOME)/src/gnuwin32/MkRules

# COLAMD/Source/Makefile uses the .c.o rule for which we need
CFLAGS=$(PKG_CPPFLAGS) -O3
# other Makefiles use these explicitly
ALL_CPPFLAGS=$(PKG_CPPFLAGS)
ALL_CFLAGS=-O3
ALL_CXXFLAGS=-O3
