/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#include "parser.hpp"

#ifdef METADYNAMICS
#include "metadynamics.hpp"

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
  char  **tmp_argv;
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
    char  **tmp_argv2;
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
    char  **tmp_argv2;
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
  return gather_runtime_errors(interp, err);
}

#endif
