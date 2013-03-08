/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2008,2009,2010 
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
/** \file adresso.c
    This is the place for adaptive resolution scheme
    Implementation of adresso.h
*/

#include "adresso.h"
#include "interaction_data.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"

/** \name Privat Functions */
/************************************************************/
/*@{*/
#ifdef ADRESS
/** calc weighting function of a distance
    @param dist distance
    @return weight of the distance
*/
double adress_wf(double dist);

#endif

/*@}*/

double adress_vars[7]       = {0, 0, 0, 0, 0, 0, 0};

#ifdef ADRESS

double adress_wf_vector(double x[3]){
  int topo=(int)adress_vars[0];
  double dist;
  int dim;
  
  
  
  switch (topo) {
  case 0:
    return 0.0;
    break;
  case 1:
    return adress_vars[1];
    break;
  case 2:
    dim=(int)adress_vars[3];
    //dist=fabs(x[dim]-adress_vars[4]);
    dist = x[dim]-adress_vars[4];
    if(dist>0)
      while(dist>box_l[dim]/2.0)
	dist = dist - box_l[dim];
    else if(dist < 0)
      while(dist< -box_l[dim]/2.0)
	dist = dist + box_l[dim];
    dist = fabs(dist);
    return adress_wf(dist);
    break;
  case 3:
    //int img_box[3];
    //double temp_pos[3];
    //for(dim=0;dim<3;dim++){
    //  img_box[dim]=0;
    //  temp_pos[dim]=x[dim];
    //}
    //fold_position(temp_pos,img_box);
    dist=distance(x,&(adress_vars[3]));
    return adress_wf(dist);
    break;
  default:
    return 0.0;
    break;
  }
}

double adress_wf(double dist){
   int wf;
   double tmp;
   
   //explicit region
   if (dist < adress_vars[1]) return 1;
   //cg regime
   else if (dist> adress_vars[1]+adress_vars[2]) return 0;
   else {
      wf=(int)adress_vars[6];
      if (wf == 0){ //cos
         tmp=PI/2/adress_vars[2]*(dist-adress_vars[1]);
         return cos(tmp)*cos(tmp);
      }
      else{ //wf == 1
         tmp=(dist-adress_vars[1]);
         return 1+2*tmp*tmp-3*tmp*tmp*tmp;
      }
   }
}

void adress_update_weights(){
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         p[i].p.adress_weight=adress_wf_vector((&p[i])->r.p);
	 //printf("LOCAL %f %f\n", p[i].r.p[0], p[i].p.adress_weight);
      }
    }
  }
  for (c = 0; c < local_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      if (ifParticleIsVirtual(&p[i])) {
         p[i].p.adress_weight=adress_wf_vector((&p[i])->r.p);
	 //printf("GHOST %f %f\n", p[i].r.p[0], p[i].p.adress_weight);
      }
    }
  }
}

#ifdef INTERFACE_CORRECTION
int adress_tab_set_params(int part_type_a, int part_type_b, char* filename)
{
    IA_parameters *data;
    FILE* fp;
    //int ic_points;
    int npoints;
    double minval,minval2, maxval, maxval2;
    int i, j, newsize;
    int token;
    double dummr;
    token = 0;
    
    data = get_ia_param_safe(part_type_a, part_type_b);
    
    if (!data)
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
    //fscanf( fp , "%d ", &ic_points);
    fscanf( fp , "%d ", &npoints);
    fscanf( fp, "%lf ", &minval);
    fscanf( fp, "%lf ", &maxval);
    
    // Set the newsize to the same as old size : only changed if a new force table is being added.
    newsize = adress_tab_forces.max;
    
    if ( data->ADRESS_TAB_npoints == 0){
        // A new potential will be added so set the number of points, the startindex and newsize
        //considering that if ic_points = 0, we have two forces: ex and cg 
        //we keep the same for npoints
        data->ADRESS_TAB_npoints    = npoints;
        data->ADRESS_TAB_startindex = adress_tab_forces.max;
        newsize += 2*npoints;
    } else {
        // We have existing data for this pair of monomer types check array sizing
        if ( data->ADRESS_TAB_npoints != npoints ) {
            fclose(fp);
            return 5;
        }
    }
    
    /* Update parameters */
    data->ADRESS_TAB_maxval    = maxval;
    data->ADRESS_TAB_minval    = minval;
    strcpy(data->ADRESS_TAB_filename,filename);
    
    /* Calculate dependent parameters */
    maxval2 = maxval*maxval;
    minval2 = minval*minval;
    data->ADRESS_TAB_stepsize = (maxval-minval)/(double)(data->ADRESS_TAB_npoints - 1);
    
    /* Allocate space for new data */
    realloc_doublelist(&adress_tab_forces,newsize);
    realloc_doublelist(&adress_tab_energies,newsize);
    
    /* Read in the new force and energy table data */
    for (i =0 ; i < npoints ; i++)
    {
        fscanf(fp,"%lf",&dummr);
        //for (j =0 ; j < ic_points + 2; j++)
        for (j =0 ; j < 2; j++)
        {
            //j = 0 -> CG FORCE
            //j = 1 -> CG_ic FORCE
            
            fscanf(fp,"%lf", &(adress_tab_forces.e[j*npoints+i+data->ADRESS_TAB_startindex]));
            fscanf(fp,"%lf", &(adress_tab_energies.e[j*npoints+i+data->ADRESS_TAB_startindex]));
        }
    }
    fclose(fp);
    
    /* broadcast interaction parameters including force and energy tables*/
    mpi_bcast_ia_params(part_type_a, part_type_b);
    
    //no force cap for the moment!
    //if (force_cap != -1.0) {
    //  mpi_cap_forces(force_cap);}
    return 0;
}

#endif

int tf_set_params(int part_type, double prefactor, char * filename)
{
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
  if (fscanf( fp , "%d ", &npoints) != 1 ||
      fscanf( fp, "%lf ", &minval) != 1 ||
      fscanf( fp, "%lf ", &maxval) != 1)
    return 5;
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
    if (fscanf(fp, "%lf", &dummr) != 1 ||
	fscanf(fp, "%lf", &(thermodynamic_forces.e[i+data->TF_TAB_startindex])) != 1 ||
	fscanf(fp, "%lf", &(thermodynamic_f_energies.e[i+data->TF_TAB_startindex])) != 1)
      return 5;
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
  
  return ES_OK;
}

#endif
