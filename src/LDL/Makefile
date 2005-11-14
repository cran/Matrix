include $(MkInclude)

OBJS = ldl.o

LIB = ../LDL.a

lib: $(LIB)

$(LIB): $(OBJS)
	$(AR) -rucs $(LIB) $(OBJS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 

ldl.o: ldl.h
