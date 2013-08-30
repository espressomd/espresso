/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
#ifndef _ADRESSO_H
#define _ADRESSO_H
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

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "virtual_sites.hpp"
#include "grid.hpp"

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

#ifdef ADRESS

/** Calc adress weight function of a vector
    @param x input vector
    @return weight of the vector
*/
double adress_wf_vector(double x[3]);

/** Calc adress weight function of a particle
    @param p input particle
    @return weight of the particle
*/
inline double adress_wf_particle(Particle *p){
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

inline double adress_non_bonded_force_weight(Particle *p1,Particle *p2){
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

inline double adress_bonded_force_weight(Particle *p1,Particle *p2){
   double weight=1.0;
  if((get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p2))->p.identity )
    weight=1.0;
  else {
    double weight_1, weight_2;
    int virtual_1, virtual_2, n_part_int=2;
    weight_1=adress_wf_particle(p1);
    virtual_1=ifParticleIsVirtual(p1);
    if( (weight_1<ROUND_ERROR_PREC) && (virtual_1==0) ) return 0;
    weight_2=adress_wf_particle(p2);
    virtual_2=ifParticleIsVirtual(p2);
    if( (weight_2<ROUND_ERROR_PREC) && (virtual_2==0) ) return 0;
    
    int sum_virtual=virtual_1+virtual_2;
    if (sum_virtual>0 && sum_virtual<n_part_int) return 0.0;
    
    weight = pow(weight_1*weight_2,2.0/((double) n_part_int));
    if(sum_virtual==n_part_int) weight=1.0-weight;
  }
  return weight;
}

inline double adress_angle_force_weight(Particle *p1, Particle *p2, Particle *p3){
  double weight=1.0;
  if((get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p2))->p.identity &&(get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p3))->p.identity )
    weight=1.0;
  else {
    double weight_1, weight_2, weight_3;
    int virtual_1, virtual_2, virtual_3, n_part_int=3;
    weight_1=adress_wf_particle(p1);
    virtual_1=ifParticleIsVirtual(p1);
    if( (weight_1<ROUND_ERROR_PREC) && (virtual_1==0) ) return 0;
    weight_2=adress_wf_particle(p2);
    virtual_2=ifParticleIsVirtual(p2);
    if( (weight_2<ROUND_ERROR_PREC) && (virtual_2==0) ) return 0;
    weight_3=adress_wf_particle(p3);
    virtual_3=ifParticleIsVirtual(p3);
    if( (weight_3<ROUND_ERROR_PREC) && (virtual_3==0) ) return 0;
    
    int sum_virtual=virtual_1+virtual_2+virtual_3;
    if (sum_virtual>0 && sum_virtual<n_part_int) return 0.0;
    
    weight = pow(weight_1*weight_2*weight_3,2.0/((double) n_part_int));
    if(sum_virtual==n_part_int) weight=1.0-weight;
  }
  return weight;
}

inline double adress_dihedral_force_weight(Particle *p1, Particle *p2, Particle *p3, Particle *p4){
  double weight=1.0;
  if((get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p2))->p.identity &&(get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p3))->p.identity && (get_mol_com_particle(p1))->p.identity == (get_mol_com_particle(p4))->p.identity)
    weight=1.0;
  else {
    double weight_1, weight_2, weight_3, weight_4;
    int virtual_1, virtual_2, virtual_3, virtual_4, n_part_int=4;
    weight_1=adress_wf_particle(p1);
    virtual_1=ifParticleIsVirtual(p1);
    if( (weight_1<ROUND_ERROR_PREC) && (virtual_1==0) ) return 0;
    weight_2=adress_wf_particle(p2);
    virtual_2=ifParticleIsVirtual(p2);
    if( (weight_2<ROUND_ERROR_PREC) && (virtual_2==0) ) return 0;
    weight_3=adress_wf_particle(p3);
    virtual_3=ifParticleIsVirtual(p3);
    if( (weight_3<ROUND_ERROR_PREC) && (virtual_3==0) ) return 0;
    weight_4=adress_wf_particle(p4);
    virtual_4=ifParticleIsVirtual(p4);
    if( (weight_4<ROUND_ERROR_PREC) && (virtual_4==0) ) return 0;
    
    int sum_virtual=virtual_1+virtual_2+virtual_3+virtual_4;
    if (sum_virtual>0 && sum_virtual<n_part_int) return 0.0;
    
    weight = pow(weight_1*weight_2*weight_3*weight_4,2.0/((double) n_part_int));
    if(sum_virtual==n_part_int) weight=1.0-weight;
  }
  return weight;
}

#ifdef INTERFACE_CORRECTION
/** Adress scheme for tabulated forces 
    - useful for particles that DO NOT 
    change their number of degrees of freedom
    and for interface pressure correction as well-
*/

inline void adress_interpolation( IA_parameters *ia_params,
				    double d[3], double dist, double force[3], int index){
  int tablepos, table_start,j;
  int inter_index = index*ia_params->ADRESS_TAB_npoints;
  double phi, dindex, fac;
  double maxval = ia_params->ADRESS_TAB_maxval;
  double minval = ia_params->ADRESS_TAB_minval;
  //int ic_points = ia_params->ADRESS_IC_npoints;
  //int max_index = 1;
  
  fac = 0.0;
  
  //if(index == max_index)
  //return;
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

inline double correction_function(double x){
  /* correction function goes between zero and one */
  double ic_s;
  ic_s = 4.0*(sqrt(x)-0.5)*(sqrt(x)-0.5);
  return ic_s;
}

///
int adress_tab_set_params(int part_type_a, int part_type_b, char* filename);

/** Adds force in an Adress way. Also useful for virial calculations */
inline void add_adress_tab_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
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

inline double inverse_weight(double w){
  if(adress_vars[0] == 2) {
    return 2/M_PI*asin(sqrt(w));
  } else {
    fprintf(stderr, "Thermodynamic force not implemented for this topology.\n");
    errexit();
  }
  return 0;
}

inline double adress_dw_dir(double pos[3], double dir[3]){
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

///
int tf_set_params(int part_type, double prefactor, char * filename);

inline double tf_profile(double x_com, int type, TF_parameters * tf_params){
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

inline void add_thermodynamic_force(Particle * p){
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

/* #endif */

#endif
/*@}*/
#endif
