/** \file initialize.c
    Implementation of \ref initialize.h "initialize.h"
*/
#include "initialize.h"
#include "global.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "binary_file.h"
#include "integrate.h"
#include "statistics.h"
#include "imd.h"
#include "random.h"
#include "communication.h"

int initialize(Tcl_Interp *interp)
{
  /*
    call the initialization of the modules here
  */
  init_random();
 
  /*
    installation of tcl commands
  */
  if(this_node==0){
  /* in integrate.c */
  Tcl_CreateCommand(interp, "integrate", integrate, 0, NULL);
  /* in global.c */
  Tcl_CreateCommand(interp, "setmd", setmd, 0, NULL);
  /* in interaction_data.c */
  Tcl_CreateCommand(interp, "inter", inter, 0, NULL);
  /* in particle_data.c */
  Tcl_CreateCommand(interp, "part", part, 0, NULL);
  /* in file binaryfile.c */
  Tcl_CreateCommand(interp, "writemd", writemd, 0, NULL);
  Tcl_CreateCommand(interp, "readmd", readmd, 0, NULL);
  /* in file statistics.c */
  Tcl_CreateCommand(interp, "mindist", mindist, 0, NULL);
  /* in file imd.c */
  Tcl_CreateCommand(interp, "imd", imd, 0, NULL);
  /* in file tcl_rand.c */
  Tcl_CreateCommand(interp, "tcl_rand", tcl_rand, 0, NULL);

  /* in interaction_data.c: make sure 0<->0 ia always exists */
  make_particle_type_exist(0);
  }
  return TCL_OK;
}
