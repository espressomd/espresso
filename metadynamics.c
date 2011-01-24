/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

#include "metadynamics.h"

/** \file metadynamics.h 
 *
 *  This file contains routines to perform metadynamics.  Right now, the
 *  reaction coordinate is defined between two particles (either distance 
 *  or z-projected distance). Note that these
 *  particles can be virtual sites, in order to handle molecules.
 *
 *  - set metadynamics options 
 *  - initialize bias forces and free energy profiles 
 *  - calculate reaction coordinate for each integration step 
 *  - apply bias force on particles
 */


#ifdef METADYNAMICS
/* metadynamics switch */
int    meta_switch       = META_OFF;
/** pid of particle 1 */
int    meta_pid1         =       -1;
/** pid of particle 2 */
int    meta_pid2         =       -1;
/** bias height */
double meta_bias_height  =    0.001;
/** bias width */
double meta_bias_width   =      0.5;

/** REACTION COORDINATE */
/** RC min */
double meta_xi_min       =        1;
/** RC max */
double meta_xi_max       =        0;
/** Force at boundaries */
double meta_f_bound      =       10;
/** Number of bins of RC */
int    meta_xi_num_bins  =      100;
double meta_xi_step      =        1;


/** Accumulated force array */
double *meta_acc_force   =     NULL;
/** Accumulated free energy profile */
double *meta_acc_fprofile=     NULL;


double *meta_cur_xi      =     NULL;
double meta_val_xi       =       0.;
double *meta_apply_direction = NULL;

int tclcommand_metadynamics(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  int err = TCL_OK;

  /* print metadynamics status */
  if(argc == 1) return tclcommand_metadynamics_print_status(interp);

  if ( ARG1_IS_S("set") )          {
    argc--;
    argv++;

    if (argc == 1) {
      Tcl_AppendResult(interp, "wrong # args: \n", (char *)NULL);
      return tclcommand_metadynamics_print_usage(interp, argc, argv);
    }
  }
  if ( ARG1_IS_S("off") )
    err = tclcommand_metadynamics_parse_off(interp, argc, argv);
  else if ( ARG1_IS_S("distance"))
    err = tclcommand_metadynamics_parse_distance(interp, argc, argv);
  else if ( ARG1_IS_S("relative_z"))
    err = tclcommand_metadynamics_parse_relative_z(interp, argc, argv);
  else if ( ARG1_IS_S("print_stat"))
    err = tclcommand_metadynamics_print_stat(interp, argc, argv);
  else if ( ARG1_IS_S("load_stat"))
    err = tclcommand_metadynamics_parse_load_stat(interp, argc, argv);
  else {
    Tcl_AppendResult(interp, "Unknown metadynamics command ", argv[1], "\n", (char *)NULL);
    return tclcommand_metadynamics_print_usage(interp, argc, argv);
  }
  return mpi_gather_runtime_errors(interp, err);
}


int tclcommand_metadynamics_print_status(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  /* metadynamics not initialized */
  if(meta_pid1 == -1 || meta_pid2 == -1) {
    Tcl_AppendResult(interp,"{ not initialized } ", (char *)NULL);
    return (TCL_OK);
  }

  /* metdynamics off */
  if(meta_switch == META_OFF) {
    Tcl_AppendResult(interp,"{ off } ", (char *)NULL);
    return (TCL_OK);
  }

  /* distance */
  if(meta_switch == META_DIST ) {
    sprintf(buffer,"%i", meta_pid1);
    Tcl_AppendResult(interp,"{ distance ",buffer, (char *)NULL);
    sprintf(buffer,"%i", meta_pid2);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_xi_min, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_xi_max, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_bias_height, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL); 
    Tcl_PrintDouble(interp, meta_bias_width, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_f_bound, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    sprintf(buffer,"%i", meta_xi_num_bins);
    Tcl_AppendResult(interp," ",buffer," } ", (char *)NULL);
  }

  /* relative_z */
  if(meta_switch == META_REL_Z ) {
    sprintf(buffer,"%i", meta_pid1);
    Tcl_AppendResult(interp,"{ relative_z ",buffer, (char *)NULL);
    sprintf(buffer,"%i", meta_pid2);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_xi_min, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_xi_max, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_bias_height, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL); 
    Tcl_PrintDouble(interp, meta_bias_width, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    Tcl_PrintDouble(interp, meta_f_bound, buffer);
    Tcl_AppendResult(interp," ",buffer, (char *)NULL);
    sprintf(buffer,"%i", meta_xi_num_bins);
    Tcl_AppendResult(interp," ",buffer," } ", (char *)NULL);
  }

  return (TCL_OK);
}

int tclcommand_metadynamics_print_usage(Tcl_Interp *interp, int argc, char **argv)
{
  Tcl_AppendResult(interp, "Usage of tcl-command metadynamics:\n", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], "' for status return or \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set off' to deactivate it \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set distance   <pid1> <pid2> <d_min> <d_max> <b_height> <b_width> <f_bound> <d_bins>' or \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " set relative_z <pid1> <pid2> <z_min> <z_max> <b_height> <b_width> <f_bound> <z_bins>' or \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " print_stat current_coord' to print value of current reaction coordinate \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " print_stat coord_values' to print all values the reaction coordinate can take \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " print_stat profile' to print the profile as a function of the reaction coordinate \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " print_stat force' to print the accumulated biased force as a function of the reaction coordinate \n ", (char *)NULL);
  Tcl_AppendResult(interp, "'", argv[0], " load_stat <profile_list> <force_list>' to load earlier simulation \n ", (char *)NULL);
  return (TCL_ERROR);
}

int tclcommand_metadynamics_parse_off(Tcl_Interp *interp, int argc, char **argv)
{
  /* set pids to -1 - invalidates the algorithm */
  meta_pid1 = -1;
  meta_pid2 = -1;
  /* switch metadynamics off */
  meta_switch = META_OFF;

  return (TCL_OK);
}

int tclcommand_metadynamics_parse_distance(Tcl_Interp *interp, int argc, char **argv)
{
  int    pid1, pid2, dbins;
  double dmin, dmax, bheight, bwidth, fbound;

  /* check number of arguments */
  if (argc < 8) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
                     argv[0]," ",argv[1]," <pid1> <pid2> <d_min> <d_max> <b_height> <b_width> <f_bound> <d_bins>\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_I(2, pid1) || !ARG_IS_I(3, pid2) || !ARG_IS_D(4, dmin) || !ARG_IS_D(5, dmax) || 
       !ARG_IS_D(6, bheight) || !ARG_IS_D(7, bwidth) || !ARG_IS_D(8, fbound) || !ARG_IS_I(9, dbins)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs two INTS, five DOUBLES, and one INT", (char *)NULL);
    return (TCL_ERROR);
  }

  if (pid1 < 0 || pid1 > max_seen_particle || pid2 < 0 || pid2 > max_seen_particle) {
    Tcl_AppendResult(interp, "pid1 and/or pid2 out of range", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (dmin < 0 || dmax < 0 || dmax < dmin || bheight < 0 || bwidth < 0 || fbound < 0 || dbins < 0) {
    Tcl_AppendResult(interp, "check parameters: inconcistency somewhere", (char *)NULL);
    return (TCL_ERROR);
  }
  
  free(meta_acc_force);
  free(meta_acc_fprofile);

  /* broadcast parameters */
  meta_pid1         = pid1;
  meta_pid2         = pid2;
  meta_bias_height  = bheight;
  meta_bias_width   = bwidth;
  meta_xi_min       = dmin;
  meta_xi_max       = dmax;
  meta_f_bound      = fbound;
  meta_xi_num_bins  = dbins;

  meta_switch = META_DIST;

  return (TCL_OK);
}


int tclcommand_metadynamics_parse_relative_z(Tcl_Interp *interp, int argc, char **argv)
{
  int    pid1, pid2, dbins;
  double dmin, dmax, bheight, bwidth, fbound;

  /* check number of arguments */
  if (argc < 8) {
    Tcl_AppendResult(interp, "wrong # args:  should be \n\"",
                     argv[0]," ",argv[1]," <pid1> <pid2> <z_min> <z_max> <b_height> <b_width> <f_bound> <z_bins>\"", (char *)NULL);
    return (TCL_ERROR);
  }

  /* check argument types */
  if ( !ARG_IS_I(2, pid1) || !ARG_IS_I(3, pid2) || !ARG_IS_D(4, dmin) || !ARG_IS_D(5, dmax) ||
       !ARG_IS_D(6, bheight) || !ARG_IS_D(7, bwidth) || !ARG_IS_D(8, fbound) || !ARG_IS_I(9, dbins)) {
    Tcl_AppendResult(interp, argv[0]," ",argv[1]," needs two INTS, five DOUBLES, and one INT", (char *)NULL);
    return (TCL_ERROR);
  }

  if (pid1 < 0 || pid1 > max_seen_particle || pid2 < 0 || pid2 > max_seen_particle) {
    Tcl_AppendResult(interp, "pid1 and/or pid2 out of range", (char *)NULL);
    return (TCL_ERROR);
  }

  if (dmax < dmin || bheight < 0 || bwidth < 0 || fbound < 0 || dbins < 0) {
    Tcl_AppendResult(interp, "check parameters: inconcistency somewhere", (char *)NULL);
    return (TCL_ERROR);
  }

  free(meta_acc_force);
  free(meta_acc_fprofile);

  /* broadcast parameters */
  meta_pid1         = pid1;
  meta_pid2         = pid2;
  meta_bias_height  = bheight;
  meta_bias_width   = bwidth;
  meta_xi_min       = dmin;
  meta_xi_max       = dmax;
  meta_f_bound      = fbound;
  meta_xi_num_bins  = dbins;

  meta_switch = META_REL_Z;

  return (TCL_OK);
}


int tclcommand_metadynamics_print_stat(Tcl_Interp *interp, int argc, char **argv)
{
  int j;
  char buffer[TCL_DOUBLE_SPACE];

  /* In case nothing has been initialized yet */
  if (meta_acc_fprofile == NULL) return (TCL_OK);

  argc -= 1; argv += 1;
  
  if ( ARG1_IS_S("current_coord") ) {
    /* Current value of the reaction coordinate */
    Tcl_PrintDouble(interp, meta_val_xi, buffer);
    Tcl_AppendResult(interp,"",buffer, (char *)NULL);
  } else if ( ARG1_IS_S("coord_values") ) {
    /* Possible values of the reaction coordinate */
    for (j = 0; j < meta_xi_num_bins; ++j) {
      Tcl_PrintDouble(interp, meta_xi_min+j*meta_xi_step, buffer);
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
    }
  } else if ( ARG1_IS_S("profile") ) {
    /* Values of the free energy profile */
    for (j = 0; j < meta_xi_num_bins; ++j) {
      Tcl_PrintDouble(interp, meta_acc_fprofile[j], buffer);
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
    }
  } else if ( ARG1_IS_S("force") ) {
    /* Values of the biased force */
    for (j = 0; j < meta_xi_num_bins; ++j) {
      Tcl_PrintDouble(interp, -1.*meta_acc_force[j], buffer);
      Tcl_AppendResult(interp,buffer," ", (char *)NULL);
    }
  } else {
    Tcl_AppendResult(interp, "Unknown option for 'metadynamics print_stat'", (char *)NULL);
    return (TCL_ERROR);
  }
  return (TCL_OK);
}


int tclcommand_metadynamics_parse_load_stat(Tcl_Interp *interp, int argc, char **argv){
  /* Parse free energy profile and biased force that were provided from an 
   * earlier simulation. Allows one to restart from a loaded state, and can 
   * even be used to allow multiple walkers communicating their data through TCL. */
  
  if(meta_switch == META_OFF) {
    Tcl_AppendResult(interp, "Metadynamics hasn't been initialized yet", (char *)NULL);
    return (TCL_ERROR);
  }
							       
  argc -= 1; argv += 1;
  
  // There should be
  if (argc != 3) {
    Tcl_AppendResult(interp, "Incorrect number of arguments: 'metadynamics load_stat <profile_list> <force_list>'", (char *)NULL);
    return (TCL_ERROR);
  }
							       
  // load free energy profile
  int i, tmp_argc, parse_error = 0, empty_line=0;
  const char  **tmp_argv;
  DoubleList profile, force;
  
  init_doublelist(&profile);
  Tcl_ResetResult(interp);
  Tcl_SplitList(interp, argv[1], &tmp_argc, &tmp_argv);
  realloc_doublelist(&profile, profile.n = tmp_argc);
  printf("profile.n %d, meta_xi_num_bins %d\n",profile.n,meta_xi_num_bins);
  /* Now check that the number of items parsed is equal to the number of bins */
  /* If there's one extra line, assume it's an empty line */
  if (profile.n == meta_xi_num_bins+1)
      empty_line = 1;
  else if (profile.n != meta_xi_num_bins) {
      Tcl_AppendResult(interp, "Size of profile list loaded is different than expected from number of bins", (char *)NULL);
      return (TCL_ERROR);
  }
  /* call meta_init() in case it has been loaded yet */
  meta_init();
  
  for(i = 0 ; i < tmp_argc-empty_line; i++) {
    int tmp_argc2;
    const char  **tmp_argv2;
    Tcl_SplitList(interp, tmp_argv[i], &tmp_argc2, &tmp_argv2);
    if (tmp_argc2 != 1) {
      Tcl_AppendResult(interp, "data set has to be a list of doubles", (char *) NULL);
      parse_error = 1; break;
    }
    if (Tcl_GetDouble(interp, tmp_argv2[0], &(profile.e[i])) == TCL_ERROR) { parse_error = 1; break; }
    /* Load data into meta_acc_fprofile */
    meta_acc_fprofile[i] = profile.e[i];
    
    Tcl_Free((char *)tmp_argv2);
  }
  Tcl_Free((char *)tmp_argv);
  if (parse_error) return TCL_ERROR;   
 

  // load force
  argc -= 1; argv += 1;
  init_doublelist(&force);
  Tcl_ResetResult(interp);
  Tcl_SplitList(interp, argv[1], &tmp_argc, &tmp_argv);
  realloc_doublelist(&force, force.n = tmp_argc);
  /* Now check that the number of items parsed is equal to the number of bins */
  if (profile.n == meta_xi_num_bins+1)
      empty_line = 1;
  else if (profile.n != meta_xi_num_bins) {
    Tcl_AppendResult(interp, "Size of force list loaded is different than expected from number of bins", (char *)NULL);
    return (TCL_ERROR);
  }
  for(i = 0 ; i < tmp_argc-empty_line; i++) {
    int tmp_argc2;
    const char  **tmp_argv2;
    Tcl_SplitList(interp, tmp_argv[i], &tmp_argc2, &tmp_argv2);
    if (tmp_argc2 != 1) {
      Tcl_AppendResult(interp, "data set has to be a list of doubles", (char *) NULL);
      parse_error = 1; break;
    }
    if (Tcl_GetDouble(interp, tmp_argv2[0], &(force.e[i])) == TCL_ERROR) { parse_error = 1; break; }
    /* Load data into meta_acc_fprofile */
    meta_acc_force[i] = -1.*force.e[i];
    
    Tcl_Free((char *)tmp_argv2);
  }
  Tcl_Free((char *)tmp_argv);
  if (parse_error) return TCL_ERROR;   

  return (TCL_OK);
}


void meta_init(){
  if(meta_switch == META_OFF)   return;


  /* Initialize arrays if they're empty. These get freed upon calling the Tcl
   * parser */
  if (meta_acc_force == NULL || meta_acc_fprofile == NULL) {
    meta_acc_force       = calloc(meta_xi_num_bins * sizeof *meta_acc_force, sizeof *meta_acc_force);
    meta_acc_fprofile    = calloc(meta_xi_num_bins * sizeof *meta_acc_fprofile, sizeof *meta_acc_fprofile);
    meta_cur_xi          = calloc(3 * sizeof *meta_cur_xi, sizeof *meta_cur_xi);
    meta_apply_direction = calloc(3 * sizeof *meta_apply_direction, sizeof *meta_apply_direction);   
  }

  /* Check that the simulation uses onle a single processor. Otherwise exit. 
  *  MPI interface *not* implemented. */
  if (n_nodes != 1) {
    char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
    ERROR_SPRINTF(errtxt,"Can't use metadynamics on more than one processor.\n");
    return;
  }

  meta_xi_step = (meta_xi_max-meta_xi_min)/(1.*meta_xi_num_bins);
}


/** Metadynamics main function:
 * - Calculate reaction coordinate
 * - Update profile and biased force
 * - apply external force
 */
void meta_perform()
{
  if(meta_switch == META_OFF)  return;

  double ppos1[3], ppos2[3], factor;
  int img1[3], img2[3], c, i, np, flag1 = 0, flag2 = 0;
  Particle *p, *p1 = NULL, *p2 = NULL;
  Cell *cell;
  
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (p[i].p.identity == meta_pid1) {
	flag1 = 1;
	p1 = &p[i];
	memcpy(ppos1, p[i].r.p, 3*sizeof(double));
	memcpy(img1, p[i].l.i, 3*sizeof(int));
	unfold_position(ppos1, img1);
	if (flag1 && flag2) {
	  /* vector r2-r1 - Not a minimal image! Unfolded position */
	  vector_subt(meta_cur_xi,ppos2,ppos1);
	  break;
	}
      }
      if (p[i].p.identity == meta_pid2) {
	flag2 = 1;
	p2 = &p[i];
	memcpy(ppos2, p[i].r.p, 3*sizeof(double));
	memcpy(img2, p[i].l.i, 3*sizeof(int));
	unfold_position(ppos2, img2);
	if (flag1 && flag2) {
	  /* vector r2-r1 - Not a minimal image! Unfolded position */
	  vector_subt(meta_cur_xi,ppos2,ppos1);
	  break;
	}
      }
    }
  }

  if (flag1 == 0 || flag2 == 0) {
      char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"Metadynamics: can't find pid1 or pid2.\n");
      return;
  }

  /* Now update free energy profile 
   * Here, we're following the functional form of 
   * Marsili etal., J Comp. Chem, 31 (2009).
   * Instead of gaussians, we use so-called Lucy's functions */

  for (i = 0; i < meta_xi_num_bins; ++i) {
    if (meta_switch == META_DIST) {
      // reaction coordinate value
      meta_val_xi = sqrt(sqrlen(meta_cur_xi));
      // Update free energy profile and biased force
      meta_acc_fprofile[i] -= calculate_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
      meta_acc_force[i] -= calculate_deriv_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);

      // direction of the bias force
      unit_vector(meta_cur_xi,meta_apply_direction);
    } else if (meta_switch == META_REL_Z) {
      // reaction coordinate value: relative height of z_pid1 with respect to z_pid2
      meta_val_xi = -1.*meta_cur_xi[2];
      // Update free energy profile and biased force
      meta_acc_fprofile[i] -= calculate_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);
      meta_acc_force[i] -= calculate_deriv_lucy(meta_xi_min+i*meta_xi_step,meta_val_xi);

      // direction of the bias force (-1 to be consistent with META_DIST: from 1 to 2)
      meta_apply_direction[0] = meta_apply_direction[1] = 0.;
      meta_apply_direction[2] = -1.;
    } else {
      char *errtxt = runtime_error(128 + 3*TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt,"Undefined metadynamics scheme.\n");
      return;      
    }
  }

  /** Apply force */

  // Calculate the strength of the applied force
  if (meta_val_xi < meta_xi_min) {
    // below the lower bound
    factor = -1. * meta_f_bound * (meta_xi_min-meta_val_xi)/meta_xi_step;
  } else if (meta_val_xi > meta_xi_max) {
    // above the upper bound
    factor = meta_f_bound * (meta_val_xi-meta_xi_max)/meta_xi_step;
  } else {
    // within the RC interval
    i = (int) dround((meta_val_xi-meta_xi_min)/meta_xi_step);
    if (i < 0) i = 0;
    if (i >= meta_xi_num_bins) i=meta_xi_num_bins-1;
    factor = meta_acc_force[i];
  }

  /* cancel previous force to external force of particle */
  for (i = 0; i<3; ++i) {
    p1->f.f[i] +=       factor * meta_apply_direction[i];
    p2->f.f[i] += -1. * factor * meta_apply_direction[i];
  }
}


/** Calculate Lucy's function */
double calculate_lucy(double xi, double xi_0)
{  
  double dist = fabs(xi-xi_0);
  if (dist <= meta_bias_width) {
    return meta_bias_height 
      * (1+2*dist/meta_bias_width) 
      * pow(1-dist/meta_bias_width,2);
  } else 
    return 0.;
}

/** Calculate derivative of Lucy function */
double calculate_deriv_lucy(double xi, double xi_0)
{
  double dist = fabs(xi-xi_0);
  if (dist <= meta_bias_width) {
    double result = -2*meta_bias_height/meta_bias_width
      * (pow(1-dist/meta_bias_width,2)
	 -(1+2*dist/meta_bias_width)*(1-dist/meta_bias_width));
    if (xi < xi_0)
      result *= -1.;
    return result;
  } else
    return 0.;
}



#endif
