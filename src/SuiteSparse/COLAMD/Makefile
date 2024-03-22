sources = Source/colamd.c Source/colamd_l.c Source/colamd_version.c
objects = $(sources:.c=.o)
archive = COLAMD.a

PKG_CPPFLAGS = -I./Include -I../SuiteSparse_config

all : $(archive)

$(archive) : $(objects)
	rm -f $@
	$(AR) -cr $@ $(objects)
	$(RANLIB) $@

.c.o :
	$(CC) $(PKG_CPPFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $*.o

clean :
	@rm -f $(objects) $(archive)
