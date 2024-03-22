sources = SuiteSparse_config.c
objects = $(sources:.c=.o)
archive = SuiteSparse_config.a

all : $(archive)

$(archive) : $(objects)
	rm -f $@
	$(AR) -cr $@ $(objects)
	$(RANLIB) $@

.c.o :
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $*.o

clean :
	@rm -f $(objects) $(archive)
