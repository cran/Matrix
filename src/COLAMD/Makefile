include ${MkInclude}

OBJS = colamd.o colamd_global.o

LIB = ../COLAMD.a

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rucs $(LIB) $(OBJS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 

colamd.o: colamd.h

