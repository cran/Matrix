sources = \
	Source/amd_1.c          Source/amd_l1.c \
	Source/amd_2.c          Source/amd_l2.c \
	Source/amd_aat.c        Source/amd_l_aat.c \
	Source/amd_control.c    Source/amd_l_control.c \
	Source/amd_defaults.c   Source/amd_l_defaults.c \
	Source/amd_dump.c       Source/amd_l_dump.c \
	Source/amd_info.c       Source/amd_l_info.c \
	Source/amd_order.c      Source/amd_l_order.c \
	Source/amd_post_tree.c  Source/amd_l_post_tree.c \
	Source/amd_postorder.c  Source/amd_l_postorder.c \
	Source/amd_preprocess.c Source/amd_l_preprocess.c \
	Source/amd_valid.c      Source/amd_l_valid.c \
	Source/amd_version.c
objects = $(sources:.c=.o)
archive = AMD.a

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
