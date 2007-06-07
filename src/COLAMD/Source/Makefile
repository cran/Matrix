include $(MkInclude)

PKG_CFLAGS = -I../Include -I../../UFconfig

LIB = ../../COLAMD.a

lib: $(LIB)

$(LIB): colamd_global.o colamd.o
	$(AR) -rucs $(LIB) colamd_global.o colamd.o

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 
