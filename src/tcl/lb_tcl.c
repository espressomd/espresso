/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file lb_tcl.c
 *
 * TCL Interface for the Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 */

#include "thermostat.h"
#include "lb_tcl.h"
#include "lb.h"
#include "parser.h"

#ifdef LB_GPU
static int lbnode_parse_set(Tcl_Interp *interp, int argc, char **argv, int *ind) {
  double f[3];
  
  while (argc > 0) {
    if(ARG0_IS_S("force")){
      if (argc < 4 ||
 	  !ARG_IS_D(1, f[0]) ||
 	  !ARG_IS_D(2, f[1]) ||
 	  !ARG_IS_D(3, f[2])
	  ) {
	Tcl_AppendResult(interp, "force expects three doubles as argument", (char *)NULL);
	return TCL_ERROR;
      }
      argc -= 4;
      argv += 4;
      if (argc > 0) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp, "Error in lbnode_extforce force. You can only change one field at the same time.", (char *)NULL);
	return ES_ERROR;
      }
    }
    else {
      Tcl_AppendResult(interp, "unknown parameter \"", argv[0], "\" to set", (char *)NULL);
      return TCL_ERROR;
    }
  }

  if (lb_lbnode_set_extforce_GPU(ind, f) == ES_ERROR) {
    Tcl_AppendResult(interp, "position is not in the LB lattice", (char *)NULL);
    return TCL_ERROR;
  }

  return ES_OK;
}

/** Parser for the \ref tclcommand_lbnode_extforce_gpu command. Can be used in future to set more values like rho,u e.g.
*/
int tclcommand_lbnode_extforce_gpu(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

  int err=ES_ERROR;
  int coord[3];

  --argc; ++argv;
  
  if (argc < 3) {
    Tcl_AppendResult(interp, "too few arguments for lbnode_extforce", (char *)NULL);
    return ES_ERROR;
  }

  if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2])) {
    Tcl_AppendResult(interp, "wrong arguments for lbnode", (char *)NULL);
    return ES_ERROR;
  } 
  argc-=3; argv+=3;

  if (argc == 0 ) { 
    Tcl_AppendResult(interp, "lbnode_extforce syntax: lbnode_extforce X Y Z [ print | set ] [ F(X) | F(Y) | F(Z) ]", (char *)NULL);
    return ES_ERROR;
  }

  if (ARG0_IS_S("set")) 
    err = lbnode_parse_set(interp, argc-1, argv+1, coord);
  else {
    Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode_extforce", (char *)NULL);
    return  ES_ERROR;
  }     
  return err;
}
#endif/* LB_GPU */

#if defined (LB) || defined (LB_GPU)
/* ********************* TCL Interface part *************************************/
/* ******************************************************************************/
void lbfluid_tcl_print_usage(Tcl_Interp *interp) {
  Tcl_AppendResult(interp, "Usage of \"lbfluid\":\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbfluid [ agrid #float ] [ dens #float ] [ visc #float ] [ tau #tau ]\n", (char *)NULL);
  Tcl_AppendResult(interp, "        [ bulk_visc #float ] [ friction #float ] [ gamma_even #float ] [ gamma_odd #float ]\n", (char *)NULL);
  Tcl_AppendResult(interp, "        [ ext_force #float #float #float ]\n", (char *)NULL);
}
void lbnode_tcl_print_usage(Tcl_Interp *interp) {
  Tcl_AppendResult(interp, "lbnode syntax:\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbnode X Y Z print [ rho | u | pi | pi_neq | boundary | populations ]\n", (char *)NULL);
  Tcl_AppendResult(interp, "     or\n", (char *)NULL);
  Tcl_AppendResult(interp, "lbnode X Y Z set [ rho | u | populations ] #nofloats", (char *)NULL);
}

/** TCL Interface: The \ref lbfluid command. */
#endif
#ifdef LB
int tclcommand_lbfluid_print_interpolated_velocity(Tcl_Interp *interp, int argc, char **argv);
#endif
int tclcommand_lbfluid(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

#if defined (LB) || defined (LB_GPU)
  argc--; argv++;

/**if we have LB the LB cpu is set by default */
#ifdef LB
  if(!(lattice_switch & LATTICE_LB_GPU)) lattice_switch = lattice_switch | LATTICE_LB;
#else
  lattice_switch = lattice_switch | LATTICE_LB_GPU;
#endif

  int err = TCL_OK;
  double floatarg;
#ifdef EXTERNAL_FORCES
  double vectarg[3];
#endif

  if (argc < 1) {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else if (ARG0_IS_S("off")) {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else if (ARG0_IS_S("init")) {
    lbfluid_tcl_print_usage(interp);
    return TCL_ERROR;
  }
  else
  	while (argc > 0) {
      if (ARG0_IS_S("gpu") || ARG0_IS_S("GPU")) {
#ifdef LB_GPU
        lattice_switch = (lattice_switch &~ LATTICE_LB) | LATTICE_LB_GPU;
        argc--; argv++;
#else
        Tcl_AppendResult(interp, "LB_GPU is not compiled in!", NULL);
        return TCL_ERROR;
#endif
      }
      else if (ARG0_IS_S("cpu") || ARG0_IS_S("CPU")) {
#ifdef LB
        lattice_switch = (lattice_switch & ~LATTICE_LB_GPU) | LATTICE_LB;
        argc--; argv++;
#else
        Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
        return TCL_ERROR;
#endif
      }
      else if (ARG0_IS_S("density") || ARG0_IS_S("dens")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "dens requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "dens must be positive", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_density(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting dens", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("grid") || ARG0_IS_S("agrid")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "agrid requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "agrid must be positive", (char *)NULL);
          return TCL_ERROR;
        } else if (0) {
          // agrid is not compatible with box_l;
          // Not necessary because this is caught on the mpi level!
        } else {
          if ( lb_lbfluid_set_agrid(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting agrid", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("tau")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "tau requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "tau must be positive", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg < time_step ) {
	        Tcl_AppendResult(interp, "tau must larger than the MD time step", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_tau(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting tau", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("viscosity") || ARG0_IS_S("visc")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "visc requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "visc must be positive", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_visc(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting viscosity", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("bulk_viscosity")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "bulk_visc requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "bulk_visc must be positive", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_bulk_visc(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting bulk_viscosity", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("friction") || ARG0_IS_S("coupling")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "friction requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (floatarg <= 0) {
	        Tcl_AppendResult(interp, "friction must be positive", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_friction(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting friction", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("ext_force")) {
#ifdef EXTERNAL_FORCES
        if ( argc < 4 || !ARG_IS_D(1, vectarg[0]) || !ARG_IS_D(2, vectarg[1]) ||  !ARG_IS_D(3, vectarg[2]) ) {
	        Tcl_AppendResult(interp, "friction requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (lb_lbfluid_set_ext_force(vectarg[0], vectarg[1], vectarg[2]) == 0) {
          argc-=4; argv+=4;
        } else {
	        Tcl_AppendResult(interp, "Unknown Error setting ext_force", (char *)NULL);
          return TCL_ERROR;
        }
      #else
        Tcl_AppendResult(interp, "External Forces not compiled in!", (char *)NULL);
         return TCL_ERROR;
      #endif
      }
      else if (ARG0_IS_S("gamma_odd")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "gamma_odd requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (fabs(floatarg) >= 1) {
	        Tcl_AppendResult(interp, "gamma_odd must < 1", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_gamma_odd(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting gamma_odd", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("gamma_even")) {
        if ( argc < 2 || !ARG1_IS_D(floatarg) ) {
	        Tcl_AppendResult(interp, "gamma_even requires 1 argument", (char *)NULL);
          return TCL_ERROR;
        } else if (fabs(floatarg) >= 1) {
	        Tcl_AppendResult(interp, "gamma_even must < 1", (char *)NULL);
          return TCL_ERROR;
        } else {
          if ( lb_lbfluid_set_gamma_even(floatarg) == 0 ) {
            argc-=2; argv+=2;
          } else {
	          Tcl_AppendResult(interp, "Unknown Error setting gamma_even", (char *)NULL);
            return TCL_ERROR;
          }
        }
      }
      else if (ARG0_IS_S("print")) {
        if ( argc < 3 || (ARG1_IS_S("vtk") && argc < 4) ) {
	        Tcl_AppendResult(interp, "lbfluid print requires at least 2 arguments. Usage: lbfluid print [vtk] velocity|boundary filename", (char *)NULL);
          return TCL_ERROR;
        }
        else {
          argc--; argv++;
          if (ARG0_IS_S("vtk")) {
          	if (ARG1_IS_S("boundary")) {
				      if ( lb_lbfluid_print_vtk_boundary(argv[2]) != 0 ) {
					      Tcl_AppendResult(interp, "Unknown Error at lbfluid print vtk boundary", (char *)NULL);
				        return TCL_ERROR;
				      }
				    }
				    else if (ARG1_IS_S("velocity")) {
				      if ( lb_lbfluid_print_vtk_velocity(argv[2]) != 0 ) {
					      Tcl_AppendResult(interp, "Unknown Error at lbfluid print vtk velocity", (char *)NULL);
				        return TCL_ERROR;
				      }
				    }
				    else {
				    	return TCL_ERROR;
				    }
				    argc-=3; argv+=3;
		      }
		      else {
		      	if (ARG0_IS_S("boundary")) {
			   	  	if ( lb_lbfluid_print_boundary(argv[1]) != 0 ) {
				    	  Tcl_AppendResult(interp, "Unknown Error at lbfluid print boundary", (char *)NULL);
			      	  return TCL_ERROR;
			      	}
			    	}
			    	else if (ARG0_IS_S("velocity")) {
			      	if ( lb_lbfluid_print_velocity(argv[1]) != 0 ) {
				    	  Tcl_AppendResult(interp, "Unknown Error at lbfluid print velocity", (char *)NULL);
			      	  return TCL_ERROR;
			      	}
			      }
				    else {
				    	return TCL_ERROR;
				    }
			      argc-=2; argv+=2;
		      }
        }
      }
      else if (ARG0_IS_S("save_ascii_checkpoint")) { 
        if (argc < 2) {
          Tcl_AppendResult(interp, "usage: lbfluid save_ascii_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } else {
          return lb_lbfluid_save_checkpoint(argv[1], 0);
        }
      }  
      else if (ARG0_IS_S("save_binary_checkpoint")) { 
        if (argc < 2) {
          Tcl_AppendResult(interp, "usage: lbfluid save_binary_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } else {
          return lb_lbfluid_save_checkpoint(argv[1], 1);
        }
      }  
      else if (ARG0_IS_S("load_ascii_checkpoint")) { 
        if (argc < 2) {
          Tcl_AppendResult(interp, "usage: lbfluid load_ascii_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } else {
          return lb_lbfluid_load_checkpoint(argv[1], 0);
        }
      }  
      else if (ARG0_IS_S("load_binary_checkpoint")) { 
        if (argc < 2) {
          Tcl_AppendResult(interp, "usage: lbfluid load_binary_checkpoint <filename>", (char *)NULL);
          return TCL_ERROR;
        } else {
          return lb_lbfluid_load_checkpoint(argv[1], 1);
        }
      }  
#ifdef LB
			else if (ARG0_IS_S("print_interpolated_velocity")) { //this has to come after print
				return tclcommand_lbfluid_print_interpolated_velocity(interp, argc-1, argv+1);
			}
#endif
      else {
    	  Tcl_AppendResult(interp, "unknown feature \"", argv[0],"\" of lbfluid", (char *)NULL);
    	  return TCL_ERROR ;
      }

      if ((err = gather_runtime_errors(interp, err)) != TCL_OK)
        return TCL_ERROR;
  }

  mpi_bcast_parameter(FIELD_LATTICE_SWITCH);

  /* thermo_switch is retained for backwards compatibility */
  thermo_switch = (thermo_switch | THERMO_LB);
  mpi_bcast_parameter(FIELD_THERMO_SWITCH);


  return TCL_OK;
#else /* !defined LB */
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}
/** Parser for the \ref tclcommand_lbnode command. */
int tclcommand_lbnode(ClientData data, Tcl_Interp *interp, int argc, char **argv) {

#if defined (LB) || defined (LB_GPU)
   int coord[3];
   int counter;
   int integer_return = 0;
   double double_return[19];

   char double_buffer[TCL_DOUBLE_SPACE];
   char integer_buffer[TCL_INTEGER_SPACE];

   for (counter = 0; counter < 19; counter++) 
     double_return[counter]=0;


   --argc; ++argv;
   if (lattice_switch & LATTICE_LB_GPU) {
   } else {
#ifdef LB
   if (lbfluid[0][0]==0) {
     Tcl_AppendResult(interp, "lbnode: lbfluid not correctly initialized", (char *)NULL);
     return TCL_ERROR;
   }
#endif
   }

   if (argc < 3) {
     lbnode_tcl_print_usage(interp);
     return TCL_ERROR;
   }

   if (!ARG_IS_I(0,coord[0]) || !ARG_IS_I(1,coord[1]) || !ARG_IS_I(2,coord[2])) {
     Tcl_AppendResult(interp, "Coordinates are not integer.", (char *)NULL);
     return TCL_ERROR;
   } 
   if (coord[0]<0 || coord[0]>box_l[0] || coord[1]<0 || coord[1]>box_l[1] || coord[2]<0 || coord[2]>box_l[2]) {
     Tcl_AppendResult(interp, "Coordinates is not a valid LB node index", (char *)NULL);
     return TCL_ERROR;
   } 
   argc-=3; argv+=3;

   if (ARG0_IS_S("print")) {
     argc--; argv++;
     while (argc > 0) {
       if (ARG0_IS_S("rho") || ARG0_IS_S("density")) {
         lb_lbnode_get_rho(coord, double_return);
         for (counter = 0; counter < 1; counter++) {
           Tcl_PrintDouble(interp, double_return[counter], double_buffer);
           Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
         }
         argc--; argv++;
       }
       else if (ARG0_IS_S("u") || ARG0_IS_S("v") || ARG0_IS_S("velocity")) {
         lb_lbnode_get_u(coord, double_return);
         for (counter = 0; counter < 3; counter++) {
           Tcl_PrintDouble(interp, double_return[counter], double_buffer);
           Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
         }
         argc--; argv++;
       }
       else if (ARG0_IS_S("pi") || ARG0_IS_S("pressure")) {
         lb_lbnode_get_pi(coord, double_return);
         for (counter = 0; counter < 6; counter++) {
           Tcl_PrintDouble(interp, double_return[counter], double_buffer);
           Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
         }
         argc--; argv++;
       }
       else if (ARG0_IS_S("pi_neq")) { /* this has to come after pi */
         lb_lbnode_get_pi_neq(coord, double_return);
         for (counter = 0; counter < 6; counter++) {
           Tcl_PrintDouble(interp, double_return[counter], double_buffer);
           Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
         }
         argc--; argv++;
       }
       else if (ARG0_IS_S("boundary")) {
         lb_lbnode_get_boundary(coord, &integer_return);
         sprintf(integer_buffer, "%d", integer_return);
				 Tcl_AppendResult(interp, integer_buffer, " ", (char *)NULL);
	 	 		 argc--; argv++;
       }
       else if (ARG0_IS_S("populations") || ARG0_IS_S("pop")) { 
         lb_lbnode_get_pop(coord, double_return);
         for (counter = 0; counter < 19; counter++) {
           Tcl_PrintDouble(interp, double_return[counter], double_buffer);
           Tcl_AppendResult(interp, double_buffer, " ", (char *)NULL);
         }
         argc--; argv++;
       }
       else {
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "unknown fluid data \"", argv[0], "\" requested", (char *)NULL);
         return TCL_ERROR;
       }
     }
   }
   else if (ARG0_IS_S("set")) {
       argc--; argv++;
       if (ARG0_IS_S("rho") || ARG0_IS_S("density")) {
         argc--; argv++;
         for (counter = 0; counter < 1; counter++) {
           if (!ARG0_IS_D(double_return[counter])) {
             Tcl_AppendResult(interp, "recieved not a double but \"", argv[0], "\" requested", (char *)NULL);
             return TCL_ERROR;
           }
           argc--; argv++;
         }
         if (lb_lbnode_set_rho(coord, double_return[0]) != 0) {
           Tcl_AppendResult(interp, "General Error on lbnode set rho.", (char *)NULL);
           return TCL_ERROR;
         }
       }
       else if (ARG0_IS_S("u") || ARG0_IS_S("v") || ARG0_IS_S("velocity")) {
         argc--; argv++;
         for (counter = 0; counter < 3; counter++) {
           if (!ARG0_IS_D(double_return[counter])) {
             Tcl_AppendResult(interp, "received not a double but \"", argv[0], "\" requested", (char *)NULL);
             return TCL_ERROR;
           }
           argc--; argv++;
         }
         if (lb_lbnode_set_u(coord, double_return) != 0) {
           Tcl_AppendResult(interp, "General Error on lbnode set u.", (char *)NULL);
           return TCL_ERROR;
         }
       }
       else if (ARG0_IS_S("pop") || ARG0_IS_S("populations") ) {
         argc--; argv++;
         for (counter = 0; counter < 19; counter++) {
           if (!ARG0_IS_D(double_return[counter])) {
             Tcl_AppendResult(interp, "recieved not a double but \"", argv[0], "\" requested", (char *)NULL);
             return TCL_ERROR;
           }
           argc--; argv++;
         }
         if (lb_lbnode_set_pop(coord, double_return) != 0) {
           Tcl_AppendResult(interp, "General Error on lbnode set pop.", (char *)NULL);
           return TCL_ERROR;
         }
       }
       else {
     Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode x y z set", (char *)NULL);
     return  TCL_ERROR;
   }
   } else {
     Tcl_AppendResult(interp, "unknown feature \"", argv[0], "\" of lbnode", (char *)NULL);
     return  TCL_ERROR;
   }
     
   return TCL_OK;
#else /* !defined LB */
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}

int tclcommand_lbfluid_print_interpolated_velocity(Tcl_Interp *interp, int argc, char **argv) {
#ifdef LB
  double p[3], v[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];
  if (argc!=3) {
    printf("usage: print_interpolated_velocity $x $y $z");
    return TCL_ERROR;
  }
  for (int i = 0; i < 3; i++) {
    if (!ARG_IS_D(i, p[i]))
      printf("usage: print_interpolated_velocity $x $y $z");
  }
  lb_lbfluid_get_interpolated_velocity_global(p, v);
  for (int i = 0; i < 3; i++) {
    Tcl_PrintDouble(interp, v[i], buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  }
  return TCL_OK;
#else
  Tcl_AppendResult(interp, "LB is not compiled in!", NULL);
  return TCL_ERROR;
#endif
}
