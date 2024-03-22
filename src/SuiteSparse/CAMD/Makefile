sources = \
	Source/camd_1.c          Source/camd_l1.c \
	Source/camd_2.c          Source/camd_l2.c \
	Source/camd_aat.c        Source/camd_l_aat.c \
	Source/camd_control.c    Source/camd_l_control.c \
	Source/camd_defaults.c   Source/camd_l_defaults.c \
	Source/camd_dump.c       Source/camd_l_dump.c \
	Source/camd_info.c       Source/camd_l_info.c \
	Source/camd_order.c      Source/camd_l_order.c \
	Source/camd_postorder.c  Source/camd_l_postorder.c \
	Source/camd_preprocess.c Source/camd_l_preprocess.c \
	Source/camd_valid.c      Source/camd_l_valid.c \
	Source/camd_version.c
objects = $(sources:.c=.o)
archive = CAMD.a

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
