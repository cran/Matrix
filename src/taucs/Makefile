include ${R_HOME}/etc/Makeconf

TAUCS_DGZ = taucs_ccs_base.c\
            taucs_sn_llt.c

TAUCS_DZ = \
  taucs_ccs_factor_llt.c\
  taucs_ccs_solve_llt.c\
  taucs_ccs_xxt.c\
  taucs_ccs_ops.c\
  taucs_vec_base.c

OBJECTS = $(TAUCS_DGZ:.c=_D.o) $(TAUCS_DGZ:.c=.o) $(TAUCS_DZ:.c=_D.o) 

TGFLAGS = -DTAUCS_CORE_GENERAL -I. 
TDFLAGS = -DTAUCS_CORE_DOUBLE -I.

LIB = ../taucs.a

lib: $(OBJECTS)
	$(AR) -rucs $(LIB) $(OBJECTS)

mostlyclean: clean

clean:
	@-rm -rf .libs _libs $(LIB)
	@-rm -f *.o 

taucs.h: taucs_private.h taucs_config_build.h taucs_config_tests.h

$(TAUCS_DZ): taucs.h
$(TAUCS_DGZ): taucs.h

$(TAUCS_DZ:.c=_D.o): %_D.o: %.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(TDFLAGS) -c $< -o $@

$(TAUCS_DGZ:.c=_D.o): %_D.o: %.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(TDFLAGS) -c $< -o $@

$(TAUCS_DGZ:.c=.o): %.o: %.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(TGFLAGS) -c $< -o $@

$(TAUCS_G:.c=.o): %.o: %.c
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(TGFLAGS) -c $< -o $@
