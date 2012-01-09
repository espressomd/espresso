/*
  Copyright (C) 2010,2011 The ESPResSo project
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

#include "tcl_interface/adresso_tcl.h"
#include "communication.h"
#include "parser.h"
#include "cells.h"

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
#endif
