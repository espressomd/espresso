// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.

#ifndef ADRESSO_H
#define ADRESSO_H
/** \file adresso.h
    This is the place for adaptive resolution scheme (adress)
    Implementation of adresso.h
    <b>Responsible:</b>
    <a href="mailto:junghans@mpip-mainz.mpg.de">Axel</a>
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

/* The list for storing the interpolation function of interface correction */
extern DoubleList ic_correction;

/** \name Exported Functions */
/************************************************************/
/*@{*/
/** Implements the Tcl command \ref tcl_adress. This allows for seetings for adress
*/
int adress_tcl(ClientData data, Tcl_Interp *interp, int argc, char **argv);

/** For the setup of the correction function, s[x] of the interface correction */
int ic(ClientData _data, Tcl_Interp *interp, int argc, char **argv);
int ic_parse(Tcl_Interp * interp, int argc, char ** argv);
int ic_read_params(char * filename);

/** #ifdef THERMODYNAMIC_FORCE */
int tf_parse(Tcl_Interp * interp, int type, double prefactor, int argc, char ** argv);
int tf_tcl(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);
/** #endif */

int manual_update_weights(ClientData _data, Tcl_Interp * interp, int argc, char ** argv);

#ifdef ADRESS
/** Calc adress weight function of a vector
    @param x[3] vector
    @return weight of the vector
*/
double adress_wf_vector(double x[3]);


/** Calc adress weight function of a particle
    @param x[3] vector
    @return weight of the particle
*/
MDINLINE double adress_wf_particle(Particle *p){
  if (p==NULL) return 0.0;
  if (ifParticleIsVirtual(p)){
    return p->p.adress_weight;
   }
  else{
    return adress_wf_particle(get_mol_com_particle(p));
   }
}
 
/** Update adress weight of all particles
*/
void adress_update_weights();

MDINLINE double adress_non_bonded_force_weight(Particle *p1,Particle *p2){
  double adress_weight_1,adress_weight_2,force_weight;
  int virtual_1,virtual_2;
  
  //NOTE this is in order of probability to appear
  adress_weight_1=adress_wf_particle(p1);
  virtual_1=ifParticleIsVirtual(p1);
  
  //if particles 1 is ex, but in the cg regime
  if ( (adress_weight_1<ROUND_ERROR_PREC) && (virtual_1==0) ) return 0.0;
  
  adress_weight_2=adress_wf_particle(p2);
  virtual_2=ifParticleIsVirtual(p2);

  //if particles 2 is ex, but in the cg regime
  if ( (adress_weight_2<ROUND_ERROR_PREC) && (virtual_2==0) ) return 0.0;

  //mixed case is captured by cg-cg interation
  if ((virtual_1+virtual_2)==1) return 0.0;

  force_weight=adress_weight_1*adress_weight_2;
  
  //both are cg
  if ((virtual_1+virtual_2)==2) {
     //both are in ex regime
     if (force_weight>1-ROUND_ERROR_PREC) return 0.0;
     force_weight=1-force_weight;
     
  }
  //both are ex -> force_weight is already set
  //if ((virtual_1+virtual_2)==0) force_weight=force_weight;
  //(ifParticleIsVirtual(p1) ==0 || ifParticleIsVirtual(p2) ==0)
  //  printf(" particle %d %d  virtual %d %d  weights  %f  %f  product %f\n", p1->p.identity, p2->p.identity, ifParticleIsVirtual(p1), ifParticleIsVirtual(p2), adress_weight_1, adress_weight_2, force_weight); 
  return force_weight;
}

MDINLINE double adress_bonded_force_weight(Particle *p1){
  //for the moment, the bonds are kept and integrated everywhere
  return 1;
  //double adress_weight_1=adress_wf_particle(p1);
  //double force_weight;
  //NOTE only ex particles have bonded interations and bonded interactions are only inside a molecule

  //particle is cg
  //if (adress_weight_1<ROUND_ERROR_PREC) return 0.0;

  //both
  //force_weight=adress_weight_1*adress_weight_1;
  //return force_weight;
}

/** Adress scheme for tabulated forces 
    - useful for particles that DO NOT 
    change their number of degrees of freedom
    and for interface pressure correction as well-
*/

MDINLINE void adress_interpolation( IA_parameters *ia_params,
				    double d[3], double dist, double force[3], int index){
  int tablepos, table_start,j;
  int inter_index = index*ia_params->ADRESS_TAB_npoints;
  double phi, dindex, fac;
  double maxval = ia_params->ADRESS_TAB_maxval;
  double minval = ia_params->ADRESS_TAB_minval;
  int ic_points = ia_params->ADRESS_IC_npoints;
  int max_index = ic_points+1;
  
  fac = 0.0;
  
  if(index == max_index)
    return;
  if ( maxval > 0 ) {
    if ( dist < maxval){ 
      table_start = ia_params->ADRESS_TAB_startindex;
      dindex = (dist-minval)/ia_params->ADRESS_TAB_stepsize;
      tablepos = (int)(floor(dindex));  
      
      if ( dist > minval ) {
       phi = dindex - tablepos;	  
       fac = adress_tab_forces.e[inter_index + table_start + tablepos]*(1-phi) + adress_tab_forces.e[inter_index + table_start + tablepos+1]*phi;
      }
      else {
	/* Use an extrapolation beyond the table */
	if ( dist > 0 ) {
	  tablepos = 0;
	  phi = dindex - tablepos;	  
	  fac = (adress_tab_forces.e[inter_index + table_start]*minval)*(1-phi) + 
	    (adress_tab_forces.e[inter_index + table_start+1]*(minval+ia_params->ADRESS_TAB_stepsize))*phi;
	  fac = fac/dist;
	}
	else { /* Particles on top of each other .. leave fac as 0.0 */
	}
      }
      
    }
    for(j=0;j<3;j++)
      force[j] += fac * d[j];
  }
}

MDINLINE double correction_function(double x, int ic_points){
  /* correction function goes between zero and one */
  double ic_s,step,dindex,phi;
  int tablepos;
  if(ic_points == 0)
    ic_s = 1;
  else {
    double ic_step = 1.0/(double)(ic_points+1); 
    double x_com = x;
    
    /* first we classify x in the corresponding bins */
    while(x_com > 0){
      x_com = x_com-ic_step;
    }
    /* now the argument x is between 0 and 1 */
    x_com = (x_com+ic_step)/ic_step;
    
    /* and we interpolate using the table , that includes the extreme points */
    step = 1.0/((double)ic_correction.max - 1.0);
    dindex = x_com/step;
    tablepos = (int)(floor(dindex));
    
    phi = dindex - tablepos;
    ic_s = ic_correction.e[tablepos]*(1.0-phi) + ic_correction.e[tablepos+1]*phi;
  }
  
  return ic_s;
}


MDINLINE int adress_tab_set_params(int part_type_a, int part_type_b, char* filename)
{
  IA_parameters *data, *data_sym;
  FILE* fp;
  int ic_points;
  int npoints;
  double minval,minval2, maxval, maxval2;
  int i, j, newsize;
  int token;
  double dummr;
  token = 0;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym)
    return 1;
  
  if (strlen(filename) > MAXLENGTH_ADRESSTABFILE_NAME-1 )
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
  fscanf( fp , "%d ", &ic_points);
  fscanf( fp , "%d ", &npoints);
  fscanf( fp, "%lf ", &minval);
  fscanf( fp, "%lf ", &maxval);
  
  // Set the newsize to the same as old size : only changed if a new force table is being added.
  newsize = adress_tab_forces.max;
  
  if ( data->ADRESS_TAB_npoints == 0){
    // A new potential will be added so set the number of points, the startindex and newsize
    //considering that if ic_points = 0, we have two forces: ex and cg 
    //we keep the same for npoints
    data->ADRESS_TAB_npoints    = data_sym->ADRESS_TAB_npoints    = npoints;
    data->ADRESS_TAB_startindex = data_sym->ADRESS_TAB_startindex = adress_tab_forces.max;
    newsize += (ic_points+2)*npoints;
  } else {
    // We have existing data for this pair of monomer types check array sizing
    if ( data->ADRESS_TAB_npoints != npoints || data->ADRESS_IC_npoints != ic_points ) {
      fclose(fp);
      return 5;
    }
  }

  /* Update parameters symmetrically */
  data->ADRESS_TAB_maxval    = data_sym->ADRESS_TAB_maxval    = maxval;
  data->ADRESS_TAB_minval    = data_sym->ADRESS_TAB_minval    = minval;
  strcpy(data->ADRESS_TAB_filename,filename);
  strcpy(data_sym->ADRESS_TAB_filename,filename);

  /* Calculate dependent parameters */
  maxval2 = maxval*maxval;
  minval2 = minval*minval;
  data->ADRESS_TAB_maxval2 = data_sym->ADRESS_TAB_maxval2 = maxval2;
  data->ADRESS_TAB_minval2 = data_sym->ADRESS_TAB_minval2 = minval2;
  data->ADRESS_TAB_stepsize = data_sym->ADRESS_TAB_stepsize = (maxval-minval)/(double)(data->ADRESS_TAB_npoints - 1);


  /* Allocate space for new data */
  realloc_doublelist(&adress_tab_forces,newsize);
  realloc_doublelist(&adress_tab_energies,newsize);

  /* Read in the new force and energy table data */
  for (i =0 ; i < npoints ; i++)
    {
      fscanf(fp,"%lf",&dummr);
      for (j =0 ; j < ic_points + 2; j++)
	{
	  //j = 0 -> CG FORCE
	  //j = ic_points + 1 -> EX FORCE
	  //the rest corresponds to the ic points
	  //that are considered EQUIDISTANT, from 
	  // w = 0 to w = 1
	  fscanf(fp,"%lf", &(adress_tab_forces.e[j*npoints+i+data->ADRESS_TAB_startindex]));
	  fscanf(fp,"%lf", &(adress_tab_energies.e[j*npoints+i+data->ADRESS_TAB_startindex]));
	}
    }
  fclose(fp);
  
  /* broadcast interaction parameters including force and energy tables*/
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  //no force cap for the moment!
  //if (tab_force_cap != -1.0) {
  //  mpi_tab_cap_forces(tab_force_cap);}
  return 0;
}

MDINLINE int adress_tab_parser(Tcl_Interp * interp,
			int part_type_a, int part_type_b,
			int argc, char ** argv)
{
  char *filename = NULL;

  /* tabulated interactions should supply a file name for a file containing
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
  case 5:
    Tcl_AppendResult(interp, "number of data or ic points does not match the existing table", (char *)NULL);
    return 0;
    
  }
  return 2;
}

/** Adds force in an Adress way. Also useful for virial calculations */
MDINLINE void add_adress_tab_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
					double d[3], double dist, double force[3])
{
  int j;
  int ic_points = ia_params->ADRESS_IC_npoints;
  int max_index = ic_points+1;
  
  double left_force[3] = {0,0,0};
  double right_force[3] = {0,0,0};
  double ex_force[3] = {0,0,0};
  double cg_force[3] = {0,0,0};
  int left_index, right_index;
  
  //ASK FOR THE WEIGHTING FUNCTIONS!!!
  double x = p1->p.adress_weight*p2->p.adress_weight;
  //double x = local_particles[p1->p.identity]->p.adress_weight*local_particles[p2->p.identity]->p.adress_weight;
  
  //EXPLICIT CASE - just for non-virtual particles
  if(x == 1){
    adress_interpolation(ia_params, d, dist, force,ic_points+1);
    return;
  }
  //COARSE-GRAINED CASE
  else if(x == 0){
    adress_interpolation(ia_params, d, dist, force, 0);
    return;
  }
  //INTERFACE PRESSURE CORRECTION: we restrict ourselves to the switching region
  else {
    //THE EXPLICIT CONTRIBUTION - just if particles are not virtual
    adress_interpolation(ia_params, d, dist, ex_force, ic_points+1);
    for(j=0;j<3;j++)
      force[j] += x*ex_force[j];
    
    //THE COARSE-GRAINED CONTRIBUTION
    //classify the position of the particle:
    if(ic_points !=0) {
      double ic_step = 1.0/((double)ic_points + 1.0);
      double w = 0;
      while(x > w+ic_step){
	left_index++;
	w = w+ic_step;
      }
      right_index = left_index+1;
      
      if(right_index < max_index){
	adress_interpolation(ia_params,d,dist,left_force,  left_index);
	adress_interpolation(ia_params,d,dist,right_force,right_index);
	
	for(j=0;j<3;j++)
	  cg_force[j] = correction_function(x,ic_points)*left_force[j] + (1.0 - correction_function(x,ic_points))*right_force[j];
      } 
      else {
	adress_interpolation(ia_params,d,dist,cg_force,left_index);
      }
    }
    else {
      //simple weighted coarse-grained force in the switching region
      adress_interpolation(ia_params,d,dist,cg_force, 0);
    }
    for(j=0;j<3;j++){
      force[j] += (1.0 - x)*cg_force[j];
    }
      return;
  }
}

/** #ifdef THERMODYNAMIC_FORCE */

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

MDINLINE double tf_profile(double x_com, int type, TF_parameters * tf_params){
  double phi, dindex, force, pol;
  int tablepos, table_start;
  //double maxval = tf_params->TF_TAB_maxval;
  double minval = tf_params->TF_TAB_minval;
  
  //if(weight == 0 || weight == 1)
  //force = 0.0;
  //else {
    table_start = tf_params->TF_TAB_startindex;
    dindex = (x_com-minval)/tf_params->TF_TAB_stepsize;
    tablepos = (int)(floor(dindex));
    phi = dindex - tablepos;
    pol = thermodynamic_forces.e[table_start + tablepos]*(1-phi) + thermodynamic_forces.e[table_start + tablepos+1]*phi;
    
    /* THERE IS NO EXTRAPOLATION!
       the table has to start and end ALWAYS at zero and one respectively
    */
    force = pol;
    //}
  
  return force;
}

MDINLINE void add_thermodynamic_force(Particle * p){
  TF_parameters *tf_params = get_tf_param(p->p.type);
  double pref   = tf_params->TF_prefactor;
  if (pref !=0){
    double weight, width, force, sign;
    int i, type;
    double dir[3] = {0,0,0};
    
    weight = p->p.adress_weight;
    if(weight>0 && weight < 1){
      type   = p->p.type;
      width  = adress_vars[2];
      sign   = adress_dw_dir(p->r.p, dir);
      
      
      force = pref*sign*tf_profile(inverse_weight(weight), type, tf_params)/width;
      
      for(i=0;i<3;i++)
	p->f.f[i] += force*dir[i];
    }
  }
}

/** #endif */

#endif
/*@}*/
#endif
