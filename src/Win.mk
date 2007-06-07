## Note: This is included from ./Makefile.win , but also
## ----  as 'MkInclude' from several */Makefile s
include $(RHOME)/src/gnuwin32/MkRules
CFLAGS=$(PKG_CFLAGS) -O2 -D__VC__
ALL_CFLAGS=$(PKG_CFLAGS) -O2

