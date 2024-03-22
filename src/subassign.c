#include "Mdefines.h"
#include "subassign.h"

enum x_slot_kind {
	x_unknown = -2,  /* NA */
	x_pattern = -1,  /*  n */
	x_double  =  0,  /*  d */
	x_logical =  1,  /*  l */
	x_integer =  2,  /*  i */
	x_complex =  3}; /*  z */

/* x[i, j] <- value  where  x=<.[gt]CMatrix>  and  value=<.sparseVector> */
#define _n_Csp_
#include "t_subassign.c"

#define _l_Csp_
#include "t_subassign.c"

#define _i_Csp_
#include "t_subassign.c"

#define _d_Csp_
#include "t_subassign.c"

#define _z_Csp_
#include "t_subassign.c"
