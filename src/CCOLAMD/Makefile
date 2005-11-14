include $(MkInclude)

OBJS = ccolamd_global.o ccolamd.o

LIB = ../CCOLAMD.a

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rucs $(LIB) $(OBJS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 

ccolamd.o: ccolamd.h
