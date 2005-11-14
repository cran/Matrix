include $(MkInclude)

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = ../Include/amd.h amd_internal.h

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))

PKG_CFLAGS = -DDINT

LIB = ../../AMD.a

lib: $(LIB)

$(LIB): amd_global.o ${AMDI}
	$(AR) -rucs $(LIB) amd_global.o ${AMDI}

amd_global.o: amd_global.c $(INC)
amd_i_%.o: amd_%.c $(INC)
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -I../Include -c $< -o $@

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 
