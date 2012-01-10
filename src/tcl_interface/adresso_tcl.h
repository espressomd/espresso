/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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
#ifndef ADRESSO_TCL_H
#define ADRESSO_TCL_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h

    For more details about adress see:
    - M. Praprotnik, L. Delle Site and K. Kremer, JCP 123, 224106, 2005. 
    - M. Praprotnik, L. Delle Site and K. Kremer, Annu. Rev. Phys. Chem. 59, 545-571, 2008. 
    - S. Poblete, M. Praprotnik, K. Kremer and L. Delle Site, J. Chem. Phys. 132, 114101, 2010. 

    For more detail about the implementation here see:
    - C. Junghans and S. Poblete, Comp. Phys. Comm. 181, 1449, 2010.
*/

#include <tcl.h>
#include "particle_data.h"
#include "virtual_sites.h"
#include "interaction_data.h"
#include "communication.h"



/** \name Exported Variables */
/************************************************************/
/*@{*/
extern double adress_vars[7];
/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command "adress". This allows for seetings for adress
*/
int tclcommand_adress(ClientData data, Tcl_Interp *interp, int argc, char **argv);

int tclcommand_update_adress_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);

#ifdef ADRESS
// This code requires the "center of mass" implementation of virtual sites
#ifndef VIRTUAL_SITES_COM
 #error Adress requires the "center of mass"-implementation  of virtual sites. Please activate it in myconfig.h
#endif
/* #ifdef THERMODYNAMIC_FORCE */
int tclcommand_thermodynamic_force_parse_opt(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv);
int tclcommand_thermodynamic_force(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);
/* #endif */

#ifdef INTERFACE_CORRECTION
/* The list for storing the interpolation function of interface correction */
//extern DoubleList ic_correction;
/** For the setup of the correction function, s[x] of the interface correction */
//int ic(ClientData _data, Tcl_Interp *interp, int argc, char **argv);
//int ic_parse(Tcl_Interp * interp, int argc, char ** argv);
//int ic_read_params(char * filename);

#endif

#ifdef INTERFACE_CORRECTION
/** Adress scheme for tabulated forces 
    - useful for particles that DO NOT 
    change their number of degrees of freedom
    and for interface pressure correction as well-
*/


/* TODO: This function is not used anywhere. To be removed?  */
MDINLINE int adress_tab_parser(Tcl_Interp * interp,
			int part_type_a, int part_type_b,
			int argc, char ** argv)
{
  char *filename = NULL;

  /* adress_tab interactions should supply a file name for a file containing
     both force and energy profiles as well as number of points, max
     values etc.
  */
  if (argc < 2) {
    Tcl_AppendResult(interp, "tabulated potentials require a filename: "
		     "<filename>",
		     (char *) NULL);
    return 0;
  }

  /* copy tabulated parameters */
  filename = argv[1];

  switch (adress_tab_set_params(part_type_a, part_type_b, filename)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "the length of the filename must be less than 256 characters,"
		     "but is \"", filename, "\"", (char *)NULL);
    return 0;
  case 3:
    Tcl_AppendResult(interp, "cannot open \"", filename, "\"", (char *)NULL);
    return 0;
  case 4:
    Tcl_AppendResult(interp, "attempt to read file \"", filename,
		     "\" failed, could not find start the start token <#>", (char *)NULL);
    return 0;
  }
  return 2;
}

/** Adds force in an Adress way. Also useful for virial calculations */
MDINLINE void add_adress_tab_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
					double d[3], double dist, double force[3])
{
  int j;
  //int ic_points = ia_params->ADRESS_IC_npoints;
  //int max_index = 1;
  
  double left_force[3] = {0,0,0};
  double right_force[3] = {0,0,0};
  //double ex_force[3] = {0,0,0};
  double cg_force[3] = {0,0,0};
  //int left_index, right_index;
  
  //ASK FOR THE WEIGHTING FUNCTIONS!!!
  double x = p1->p.adress_weight*p2->p.adress_weight;
  //double x = local_particles[p1->p.identity]->p.adress_weight*local_particles[p2->p.identity]->p.adress_weight;
  
  //NO EXPLICIT CASE!!!
  //EXPLICIT CASE - just for non-virtual particles
  if(x == 1){
    //adress_interpolation(ia_params, d, dist, force,ic_points+1);
  return;
  }
  //COARSE-GRAINED CASE
  else if(x == 0){
    adress_interpolation(ia_params, d, dist, cg_force, 0);
    for(j=0;j<3;j++){
      force[j] += cg_force[j];
    }
    //if (sqrt(cg_force[0]*cg_force[0]+cg_force[1]*cg_force[1]+cg_force[2]*cg_force[2])!=0)
    //printf("%f   %f\n", dist, sqrt(cg_force[0]*cg_force[0]+cg_force[1]*cg_force[1]+cg_force[2]*cg_force[2]));
    return;
  }
  //INTERFACE PRESSURE CORRECTION: we restrict ourselves to the switching region
  else {
    //THE EXPLICIT CONTRIBUTION - just if particles are not virtual
    //adress_interpolation(ia_params, d, dist, ex_force, ic_points+1);
    //for(j=0;j<3;j++)
    // force[j] += x*ex_force[j];
    
    //THE COARSE-GRAINED CONTRIBUTION
    //classify the position of the particle:
    //if(ic_points !=0) {
    //double ic_step = 1.0/((double)ic_points + 1.0);
    // double w = 0;
    // while(x > w+ic_step){
    //left_index++;
    //w = w+ic_step;
    //}
    //right_index = left_index+1;
    
    //if(right_index < max_index){
    adress_interpolation(ia_params,d,dist,left_force,  0);
    adress_interpolation(ia_params,d,dist,right_force, 1);
    
    for(j=0;j<3;j++)
      cg_force[j] = correction_function(x)*left_force[j] + (1.0 - correction_function(x))*right_force[j];
    //}       else {
    //adress_interpolation(ia_params,d,dist,cg_force,left_index);
    //}
  
    for(j=0;j<3;j++){
      force[j] += cg_force[j];
    }
    return;
  }
  
}
#endif

/* #ifdef THERMODYNAMIC_FORCE */

MDINLINE double inverse_weight(double w){
  if(adress_vars[0] == 2) {
    return 2/M_PI*asin(sqrt(w));
  } else {
    fprintf(stderr, "Thermodynamic force not implemented for this topology.\n");
    errexit();
  }
  return 0;
}

MDINLINE double adress_dw_dir(double pos[3], double dir[3]){
  int topo=(int)adress_vars[0];
  double dist, mod=0;
  int i, dim;
  
  for(i=0;i<3;i++)
    dir[i]=0.0;
  
  switch (topo) {
  case 0:
    return 0.0;
    break;
  case 1:
    return 0.0;
    break;
  case 2:
    dim=(int)adress_vars[3];
    //dist=fabs(x[dim]-adress_vars[4]);
    dist = pos[dim]-adress_vars[4];
    if(dist>0)
      while(dist>box_l[dim]/2.0)
	dist = dist - box_l[dim];
    else if(dist < 0)
      while(dist< -box_l[dim]/2.0)
	dist = dist + box_l[dim];
    dir[dim]=1;
    if(dist>0)
      return -1;
    else return 1;
    
    break;
  case 3:
    /* NOT TESTED */
    dist=distance(pos,&(adress_vars[3]));
    for(i=0;i<3;i++)
      mod += (pos[i]-adress_vars[3+i])*(pos[i]-adress_vars[3+i]);
    if(mod == 0){
      fprintf(stderr,"Particle located at the center of the box: Thermodynamic force not defined.\n");
      errexit();
    }
    for(i=0;i<3;i++)
      dir[i]=(pos[i]-adress_vars[3+i])/mod;
    if(dist < adress_vars[1]+adress_vars[2])
      return -1;
    else return 1;
    break;
  default:
    return 0.0;
    break;
  }
}

MDINLINE int tf_set_params(int part_type, double prefactor, char * filename){
  TF_parameters *data;
  FILE *fp;
  int npoints;
  double minval, maxval;
  int i, newsize;
  int token = 0;
  double dummr;
  
  make_particle_type_exist(part_type);
  data = get_tf_param(part_type);
  if (!data)
    return 1;
  
  if (strlen(filename) > MAXLENGTH_TABFILE_NAME-1 )
    return 2;
  
  /*Open the file containing force and energy tables */
  fp = fopen( filename , "r");
  if ( !fp )
    return 3;
  
  /*Look for a line starting with # */
  while ( token != EOF) {
    token = fgetc(fp);
    if ( token == 35 ) { break; } // magic number for # symbol
  }
  if ( token == EOF ) {
    fclose(fp);
    return 4;
  }
  
  /* First read two important parameters we read in the data later*/
  fscanf( fp , "%d ", &npoints);
  fscanf( fp, "%lf ", &minval);
  fscanf( fp, "%lf ", &maxval);
  // Set the newsize to the same as old size : only changed if a new force table is being added.
  newsize = thermodynamic_forces.max;
  if ( data->TF_TAB_npoints == 0){
    // A new potential will be added so set the number of points, the startindex and newsize
    data->TF_TAB_npoints = npoints;
    data->TF_TAB_startindex = thermodynamic_forces.max;
    newsize += npoints;
  } else {
    // We have existing data for this pair of monomer type check array sizing
    if ( data->TF_TAB_npoints != npoints ) {
      fclose(fp);
      return 5;
    }
  }
  
  /* Update parameters */
  data->TF_TAB_maxval = maxval;
  data->TF_TAB_minval = minval;
  strcpy(data->TF_TAB_filename, filename);
  data->TF_prefactor = prefactor;
  
  data->TF_TAB_stepsize = (maxval-minval)/(double)(data->TF_TAB_npoints - 1);
  
  /* Allocate space for new data */
  realloc_doublelist(&thermodynamic_forces, newsize);
  realloc_doublelist(&thermodynamic_f_energies, newsize);
  
  /* Read in the new force and energy table data */
  for (i = 0 ; i < npoints ; i++){
    fscanf(fp, "%lf", &dummr);
    fscanf(fp, "%lf", &(thermodynamic_forces.e[i+data->TF_TAB_startindex]));
    fscanf(fp, "%lf", &(thermodynamic_f_energies.e[i+data->TF_TAB_startindex]));
    if(i==0 && dummr !=0) {
      fprintf(stderr, "First point of the thermodynamic force has to be zero.\n");
      errexit();
    }
    else if (i== npoints-1 && dummr != 1){
      fprintf(stderr, "Last point of the thermodynamic force has to be one.\n");
      errexit();
    }
  }
  
  fclose(fp);
  
  /* broadcast interaction parameters including force and energy tables */
  mpi_bcast_tf_params(part_type);
  
  return TCL_OK;
}

/* #endif */

#endif
/*@}*/
#endif
