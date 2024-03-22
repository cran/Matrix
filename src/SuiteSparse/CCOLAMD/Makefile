sources = Source/ccolamd.c Source/ccolamd_l.c Source/ccolamd_version.c
objects = $(sources:.c=.o)
archive = CCOLAMD.a

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
