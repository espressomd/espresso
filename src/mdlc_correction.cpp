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
/** \file p3m.hpp  code for calculating the MDLC (magnetic dipolar layer correction).
 *  Purpose:   get the corrections for dipolar 3D algorithms 
 *             when applied to a slab geometry and dipolar
 *	      particles. DLC & co
 *  Article:   A. Brodka, Chemical Physics Letters 400, 62-67 (2004).
 * 	      
 *	      We also include a tuning function that returns the
 *	      cut-off necessary to attend a certain accuracy.
 *	      
 *  Restrictions: the slab must be such that the z is the short 
 *                direction. Othewise we get trash.    	      
 * 
 */

#include "utils.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "domain_decomposition.hpp"
#include "particle_data.hpp"
#include "communication.hpp"
#include "p3m-dipolar.hpp"
#include "cells.hpp"
#include "mdlc_correction.hpp"

#ifdef DIPOLES
   
DLC_struct dlc_params = { 1e100, 0, 0, 0, 0};

static int n_local_particles=0;
static double mu_max;

void calc_mu_max() {
  Cell *cell;
  Particle *part;
  int i,c,np;

  mu_max = -1;

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if(mu_max <  part[i].p.dipm ) {  mu_max=part[i].p.dipm;}
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &mu_max, 1, MPI_DOUBLE, MPI_MAX, comm_cart);
}

/* ******************************************************************* */

inline double g1_DLC_dip(double g,double x) {
  double a,c,cc2,x3;
  c=g/x;
  cc2=c*c;
  x3=x*x*x;
  a=g*g*g/x+1.5*cc2+1.5*g/x3+0.75/(x3*x);
  return a;
}	
/* ******************************************************************* */


inline double g2_DLC_dip(double g,double x) {
  double a,x2;
  x2=x*x;
  a=g*g/x+2.0*g/x2+2.0/(x2*x); 
  return a;
}	
/* ******************************************************************* */    

/* Subroutine designed to  compute Mx, My, Mz and Mtotal  */
     
double slab_dip_count_mu(double *mt, double *mx, double *my)
{  
  Cell *cell;
  Particle *part;
  int i,c,np;
  double node_sums[3], tot_sums[3],Mz,M,My,Mx;
 

  node_sums[0]=0.0; tot_sums[0]=0.0;
  node_sums[1]=0.0; tot_sums[1]=0.0;
  node_sums[2]=0.0; tot_sums[2]=0.0;
 
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if( part[i].p.dipm != 0.0 ) {
	node_sums[0] +=part[i].r.dip[0];	 	 
	node_sums[1] +=part[i].r.dip[1];	 	 
	node_sums[2] +=part[i].r.dip[2];	 	 
      }
    }
  }
   
  MPI_Allreduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, comm_cart);

  M= sqrt(tot_sums[0]*tot_sums[0]+tot_sums[1]*tot_sums[1]+tot_sums[2]*tot_sums[2]);
  Mz=tot_sums[2]; 
  Mx=tot_sums[0]; 
  My=tot_sums[1]; 
   
  //fprintf(stderr,"Mz=%20.15le \n",Mz);
  //fprintf(stderr,"M=%20.15le \n",M);
 
  *mt=M;
  *mx=Mx;
  *my=My;

  return Mz;
}    
/* ******************************************************************* */    
     
     
/* ****************************************************************************************************
   Compute the dipolar DLC corrections for forces and torques.
   Algorithm implemented accordingly to the paper of A. Brodka, Chem. Phys. Lett. 400, 62-67, (2004).
   ****************************************************************************************************
   */

double get_DLC_dipolar(int kcut,double *fx, double *fy, double *fz, double *tx, double *ty, double *tz){

  int    ix,iy,ip; 
  double gx,gy,gr;

  double S[4] = {0.0,0.0,0.0,0.0}; // S of Brodka methode, oder is S[4] = {Re(S+), Im(S+), Re(S-), Im(S-)}
  double *ReSjp=NULL,*ReSjm=NULL;
  double *ImSjp=NULL,*ImSjm=NULL;
  double *ReGrad_Mup=NULL,*ImGrad_Mup=NULL;
  double *ReGrad_Mum=NULL,*ImGrad_Mum=NULL;
  double a,b,c,d,er,ez,f,fa1;
  double s1,s2,s3,s4;
  double s1z,s2z,s3z,s4z;
  double ss;
  double energy,piarea,facux,facuy;
  int cc,j,np;
  Cell *cell = NULL;   
  Particle *p1;

  facux=2.0*M_PI/box_l[0];
  facuy=2.0*M_PI/box_l[1];
      
  energy=0.0;
    
  ReSjp= (double *) malloc(sizeof(double)*n_local_particles);
  ReSjm= (double *) malloc(sizeof(double)*n_local_particles);
  ImSjp= (double *) malloc(sizeof(double)*n_local_particles);
  ImSjm= (double *) malloc(sizeof(double)*n_local_particles);
  ReGrad_Mup = (double *) malloc(sizeof(double)*n_local_particles);
  ImGrad_Mup = (double *) malloc(sizeof(double)*n_local_particles);
  ReGrad_Mum = (double *) malloc(sizeof(double)*n_local_particles);
  ImGrad_Mum = (double *) malloc(sizeof(double)*n_local_particles);

  for(ix=-kcut;ix<=+kcut;ix++){
    for(iy=-kcut;iy<=+kcut;iy++){
      if(!(ix==0 && iy==0)){  
	gx=(double)ix*facux;
	gy=(double)iy*facuy;
	  
	gr=sqrt(gx*gx+gy*gy);
	  
	fa1=1./(gr*(exp(gr*box_l[2])-1.0));   //We assume short slab direction is z direction
	  
	// ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g
	  
	S[0] = S[1] = S[2] = S[3] = 0.0;
	  
	ip=0;
	for(cc=0; cc<local_cells.n;cc++){      
	  cell = local_cells.cell[cc];
	  p1   = cell->part;
	  np   = cell->n;

	  for(j = 0; j < np; j++)  {
	    if(p1[j].p.dipm>0){
	
	      a=gx*p1[j].r.dip[0]+gy*p1[j].r.dip[1];
	      b=gr*p1[j].r.dip[2];
	      er=gx*p1[j].r.p[0] +gy*p1[j].r.p[1] ;
	      ez=gr*p1[j].r.p[2];
	      c=cos(er);
	      d=sin(er);
	      f=exp(ez);
		
	      ReSjp[ip]=(b*c-a*d)*f;
	      ImSjp[ip]=(c*a+b*d)*f;
	      ReSjm[ip]=(-b*c-a*d)/f;
	      ImSjm[ip]=(c*a-b*d)/f;
	      ReGrad_Mup[ip]=c*f;
	      ReGrad_Mum[ip]=c/f;
	      ImGrad_Mup[ip]=d*f;
	      ImGrad_Mum[ip]=d/f;
		
	      S[0]+=ReSjp[ip];
	      S[1]+=ImSjp[ip];
	      S[2]+=ReSjm[ip];
	      S[3]+=ImSjm[ip];
	    }  
	    ip++;
	  }      
	}      

	MPI_Allreduce(MPI_IN_PLACE, S, 4, MPI_DOUBLE, MPI_SUM, comm_cart);

	//We compute the contribution to the energy ............
	     
	//s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!
	          
	energy+=fa1*((S[0]*S[2]+S[1]*S[3])*2.0); 
     
	// ... Now we can compute the contributions to E,Fj,Ej for the current g-value
	ip=0;
	for(cc=0; cc<local_cells.n;cc++){      
	  cell = local_cells.cell[cc];
	  p1   = cell->part;
	  np   = cell->n;

	  for(j = 0; j < np; j++)  {
	    if(p1[j].p.dipm>0){
	      //We compute the contributions to the forces ............ 

	      s1=-(-ReSjp[ip]*S[3]+ImSjp[ip]*S[2]); 
	      s2=+( ReSjm[ip]*S[1]-ImSjm[ip]*S[0]); 
	      s3=-(-ReSjm[ip]*S[1]+ImSjm[ip]*S[0]); 
	      s4=+( ReSjp[ip]*S[3]-ImSjp[ip]*S[2]); 
	
	      s1z=+(ReSjp[ip]*S[2]+ImSjp[ip]*S[3]);
	      s2z=-(ReSjm[ip]*S[0]+ImSjm[ip]*S[1]);    
	      s3z=-(ReSjm[ip]*S[0]+ImSjm[ip]*S[1]);    
	      s4z=+(ReSjp[ip]*S[2]+ImSjp[ip]*S[3]);
	    
	      ss=s1+s2+s3+s4;
	      fx[ip]+=fa1*gx*ss;
	      fy[ip]+=fa1*gy*ss;
	    
	      fz[ip]+=fa1*gr*(s1z+s2z+s3z+s4z);
	     
	      //We compute the contributions to the electrical field ............
	  
	      s1=-(-ReGrad_Mup[ip]*S[3]+ImGrad_Mup[ip]*S[2]); 
	      s2=+( ReGrad_Mum[ip]*S[1]-ImGrad_Mum[ip]*S[0]); 
	      s3=-(-ReGrad_Mum[ip]*S[1]+ImGrad_Mum[ip]*S[0]); 
	      s4=+( ReGrad_Mup[ip]*S[3]-ImGrad_Mup[ip]*S[2]); 
	
	      s1z=+(ReGrad_Mup[ip]*S[2]+ImGrad_Mup[ip]*S[3]);
	      s2z=-(ReGrad_Mum[ip]*S[0]+ImGrad_Mum[ip]*S[1]); 
	      s3z=-(ReGrad_Mum[ip]*S[0]+ImGrad_Mum[ip]*S[1]); 
	      s4z=+(ReGrad_Mup[ip]*S[2]+ImGrad_Mup[ip]*S[3]);
	    
	      ss=s1+s2+s3+s4;
	      tx[ip]+=fa1*gx*ss;
	      ty[ip]+=fa1*gy*ss;
	    
	      tz[ip]+=fa1*gr*(s1z+s2z+s3z+s4z);
	    }//if dipm>0 ....
	    ip++;
          }//loop j  
	} //loop cc
      
      }//end of if(ii> ... 
    
    }} //end of loops for gx,gy
 
 
  //Convert from the corrections to the Electrical field to the corrections for the torques ....
  
  
  //printf("Electical field: Ex %le, Ey %le, Ez %le",-tx[0]*M_PI/(box_l[0]*box_l[1]),-ty[0]*M_PI/(box_l[0]*box_l[1]),
  //-tz[0]*M_PI/(box_l[0]*box_l[1])  );
  
  ip=0;
  for(cc=0; cc<local_cells.n;cc++){							  
    cell = local_cells.cell[cc];							  
    p1   = cell->part;								  
    np   = cell->n;									  

    for(j = 0; j < np; j++)  {
      if(p1[j].p.dipm>0){
	a=p1[j].r.dip[1]*tz[ip]-p1[j].r.dip[2]*ty[ip];
	b=p1[j].r.dip[2]*tx[ip]-p1[j].r.dip[0]*tz[ip];
	c=p1[j].r.dip[0]*ty[ip]-p1[j].r.dip[1]*tx[ip];
	tx[ip]=a;
	ty[ip]=b;
	tz[ip]=c;
	 
      }
      ip++;
    }
  }
 
  // Multiply by the factors we have left during the loops

  //printf("box_l: %le %le %le \n",box_l[0],box_l[1],box_l[2]);

  piarea=M_PI/(box_l[0]*box_l[1]);

  for(j = 0; j < n_local_particles; j++)  {
    fx[j]*=piarea;
    fy[j]*=piarea;
    fz[j]*=piarea;
    tx[j]*=piarea;
    ty[j]*=piarea;
    tz[j]*=piarea;
  }  
 
  energy*=(-piarea);
  
  //fclose(FilePtr);

  free(ReSjp);
  free(ReSjm);
  free(ImSjp);
  free(ImSjm);
  free(ReGrad_Mup);
  free(ImGrad_Mup);
  free(ReGrad_Mum);
  free(ImGrad_Mum);


  //printf("Energy0= %20.15le \n",energy);


  return energy;
}    
/* ******************************************************************* */


/* ****************************************************************************************************
   Compute the dipolar DLC corrections 
   Algorithm implemented accordingly to the paper of A. Brodka, Chem. Phys. Lett. 400, 62-67, (2004).
   ****************************************************************************************************
   */

double get_DLC_energy_dipolar(int kcut){

  int    ix,iy,ip; 
  double gx,gy,gr;

  double S[4], global_S[4];
  double a,b,c,d,er,ez,f,fa1;
  double s1;
  double energy,piarea,facux,facuy;
  int    cc,j,np;  										    
  Cell   *cell = NULL;   
  Particle *p1;  									    

  //FILE *FilePtr;
  // char File_Name[40];

  n_local_particles = 0;
  for(cc=0; cc<local_cells.n;cc++) {
    cell = local_cells.cell[cc];
    n_local_particles += cell->n;
  }

   
  facux=2.0*M_PI/box_l[0];
  facuy=2.0*M_PI/box_l[1];
  
  energy=0.0;

  for(ix=-kcut;ix<=+kcut;ix++){ 	
    for(iy=-kcut;iy<=+kcut;iy++){ 	

            
      if(!(ix==0 && iy==0)){  
	gx=(double)ix*facux;
	gy=(double)iy*facuy;
       
	gr=sqrt(gx*gx+gy*gy);
      
        fa1=1./(gr*(exp(gr*box_l[2])-1.0));   //We assume short slab direction is z direction
       
	// ... Compute S+,(S+)*,S-,(S-)*, and Spj,Smj for the current g
            
	S[0] = S[1] = S[2] = S[3] = 0.0;
      
	ip=0;
	for(cc=0; cc<local_cells.n;cc++){
	  cell = local_cells.cell[cc];
	  p1   = cell->part;
	  np   = cell->n;

	  for(j = 0; j < np; j++)  {
	    if(p1[j].p.dipm>0){
	    
	      a=gx*p1[j].r.dip[0]+gy*p1[j].r.dip[1];{
		b=gr*p1[j].r.dip[2];
		er=gx*p1[j].r.p[0] +gy*p1[j].r.p[1] ;
		ez=gr*p1[j].r.p[2];
		c=cos(er);
		d=sin(er);
		f=exp(ez);
	    
		S[0]+=(b*c-a*d)*f;
		S[1]+=(c*a+b*d)*f;
		S[2]+=(-b*c-a*d)/f;
		S[3]+=(c*a-b*d)/f;
	      }
	   
	    }
	    ip++;
	  } 										      
	}											      
	MPI_Reduce(S, global_S, 4, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
      
	//We compute the contribution to the energy ............
	s1=(global_S[0]*global_S[2]+global_S[1]*global_S[3]);
	//s2=(ReSm*ReSp+ImSm*ImSp); s2=s1!!!
	     
	energy+=fa1*(s1*2.0); 
 
      }//end of if(... 
    }} //end of loops for gx,gy
 
  
  // Multiply by the factors we have left during the loops
 
  piarea=M_PI/(box_l[0]*box_l[1]);
  energy*=(-piarea);
  return (this_node == 0) ? energy : 0.0;
}    
/* ***************************************************************** */
     

     

/* **************************************************************************
********** Compute and add the terms needed to correct the 3D dipolar*****
********** methods when we have an slab geometry *************************
************************************************************************** */
     
void    add_mdlc_force_corrections(){
  Cell *cell;
  Particle *p;
  int i,c,np,ip;
  int cc;
  int dip_DLC_kcut;
  double  *dip_DLC_f_x=NULL,*dip_DLC_f_y=NULL,*dip_DLC_f_z=NULL;
  double  *dip_DLC_t_x=NULL,*dip_DLC_t_y=NULL,*dip_DLC_t_z=NULL;
  double  dip_DLC_energy=0.0;
  double  mz=0.0,mx=0.0,my=0.0,volume,mtot=0.0;
#if defined(ROTATION) && defined(DP3M)
  double dx,dy,dz,correps;
#endif  

  dip_DLC_kcut=dlc_params.far_cut ;

  n_local_particles = 0;
  for(cc=0; cc<local_cells.n;cc++) {
    cell = local_cells.cell[cc];
    n_local_particles += cell->n;
  }
          
  volume=box_l[0]*box_l[1]*box_l[2];
	
  // --- Create arrays that should contain the corrections to
  //     the forces and torques, and set them to zero.   
 
  dip_DLC_f_x = (double *) malloc(sizeof(double)*n_part);
  dip_DLC_f_y = (double *) malloc(sizeof(double)*n_part);
  dip_DLC_f_z = (double *) malloc(sizeof(double)*n_part);
	 
  dip_DLC_t_x = (double *) malloc(sizeof(double)*n_part);
  dip_DLC_t_y = (double *) malloc(sizeof(double)*n_part);
  dip_DLC_t_z = (double *) malloc(sizeof(double)*n_part);


  for(i=0;i<n_local_particles;i++){
    dip_DLC_f_x[i] = 0.0;
    dip_DLC_f_y[i] = 0.0;
    dip_DLC_f_z[i] = 0.0;
	 
    dip_DLC_t_x[i] = 0.0;
    dip_DLC_t_y[i] = 0.0;
    dip_DLC_t_z[i] = 0.0;
  }   

  //---- Compute the corrections ----------------------------------  
     
  //First the DLC correction  
  dip_DLC_energy+=coulomb.Dprefactor*get_DLC_dipolar(dip_DLC_kcut,dip_DLC_f_x,dip_DLC_f_y,dip_DLC_f_z,dip_DLC_t_x,dip_DLC_t_y,dip_DLC_t_z);    
 
  //            printf("Energy DLC                                  = %20.15le \n",dip_DLC_energy);
 
  //Now we compute the the correction like Yeh and Klapp to take into account the fact that you are using a
  //3D PBC method which uses spherical summation instead of slab-wise sumation. Slab-wise summation is the one 
  //required to apply DLC correction.  This correction is often called SDC = Shape Dependent Correction.
  //See Brodka, Chem. Phys. Lett. 400, 62, (2004).
	   
  mz=slab_dip_count_mu(&mtot, &mx, &my);
#ifdef DP3M
  if(coulomb.Dmethod == DIPOLAR_MDLC_P3M) {
    if(dp3m.params.epsilon == P3M_EPSILON_METALLIC) {	
      dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz);
    }
    else{   
      dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz-mtot*mtot/(2.0*dp3m.params.epsilon+1.0)); 
    }
  }
  else
#endif
    {
      dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz);
      fprintf(stderr,"You are not using the P3M method, therefore p3m.epsilon is unknown, I assume metallic borders \n");   
    }  
	         
  // --- Transfer the computed corrections to the Forces, Energy and torques
  //	of the particles
	 
  ip = 0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p    = cell->part;
    np   = cell->n;
    for(i=0; i<np; i++) { 
      if( (p[i].p.dipm) != 0.0 ) {

	p[i].f.f[0] += coulomb.Dprefactor*dip_DLC_f_x[ip];    
	p[i].f.f[1] += coulomb.Dprefactor*dip_DLC_f_y[ip];
	p[i].f.f[2] += coulomb.Dprefactor*dip_DLC_f_z[ip]; //SDC correction term is zero for the forces
   
#if defined(ROTATION) && defined(DP3M)
    double correc= 4.*M_PI/volume;
	//in the Next lines: the second term (correc*...)is the SDC correction for the torques
	if(dp3m.params.epsilon == P3M_EPSILON_METALLIC) {	

	  dx=0.0;
	  dy=0.0;
	  dz=correc*(-1.0)*mz;    
		
	  p[i].f.torque[0] +=coulomb.Dprefactor*(dip_DLC_t_x[ip]+p[i].r.dip[1]*dz  - p[i].r.dip[2]*dy ) ; 
	  p[i].f.torque[1] +=coulomb.Dprefactor*(dip_DLC_t_y[ip]+p[i].r.dip[2]*dx  - p[i].r.dip[0]*dz ) ;
	  p[i].f.torque[2] +=coulomb.Dprefactor*(dip_DLC_t_z[ip]+p[i].r.dip[0]*dy  - p[i].r.dip[1]*dx ); 


	}else{
		
	  correps= correc/(2.0*dp3m.params.epsilon+1.0);
	  dx=correps*mx;
	  dy=correps*my;
	  dz=correc*(-1.0+1./(2.0*dp3m.params.epsilon+1.0))*mz;    
		
	  p[i].f.torque[0] +=coulomb.Dprefactor*(dip_DLC_t_x[ip]+p[i].r.dip[1]*dz  - p[i].r.dip[2]*dy ) ; 
	  p[i].f.torque[1] +=coulomb.Dprefactor*(dip_DLC_t_y[ip]+p[i].r.dip[2]*dx  - p[i].r.dip[0]*dz ) ;
	  p[i].f.torque[2] +=coulomb.Dprefactor*(dip_DLC_t_z[ip]+p[i].r.dip[0]*dy  - p[i].r.dip[1]*dx ); 
		
		   
	}   
#endif  	      
      }  
      ip++;
    }	 
  }
	
  //--- Free the memory used for computing the corrections ----------------    	
      
  free(dip_DLC_f_x);
  free(dip_DLC_f_y);
  free(dip_DLC_f_z);
  free(dip_DLC_t_x);
  free(dip_DLC_t_y);
  free(dip_DLC_t_z);
   
    

}    
/* ***************************************************************** */
     


/* **************************************************************************
********** Compute and add the terms needed to correct the energy of *****
********** 3D dipolar methods when we have an slab geometry          *****
************************************************************************** */
   
double     add_mdlc_energy_corrections(){
  double  dip_DLC_energy=0.0;
  double  mz=0.0,mx=0.0,my=0.0,volume,mtot=0.0;
  int dip_DLC_kcut;
     
  volume=box_l[0]*box_l[1]*box_l[2];
     
  dip_DLC_kcut=dlc_params.far_cut ;
     
  //---- Compute the corrections ----------------------------------  
     
  //First the DLC correction  
  dip_DLC_energy+=coulomb.Dprefactor*get_DLC_energy_dipolar(dip_DLC_kcut);    
     
  //           printf("Energy DLC                                  = %20.15le \n",dip_DLC_energy);
     
  //Now we compute the the correction like Yeh and Klapp to take into account the fact that you are using a
  //3D PBC method which uses spherical summation instead of slab-wise sumation. Slab-wise summation is the one 
  //required to apply DLC correction.  This correction is often called SDC = Shape Dependent Correction.
  //See Brodka, Chem. Phys. Lett. 400, 62, (2004).
     
  mz=slab_dip_count_mu(&mtot, &mx, &my);
     
  if(this_node == 0) {
#ifdef DP3M      
    if(coulomb.Dmethod == DIPOLAR_MDLC_P3M) {
      if(dp3m.params.epsilon == P3M_EPSILON_METALLIC) {
	dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz);
      }
      else{   
	dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz-mtot*mtot/(2.0*dp3m.params.epsilon+1.0)); 
      }
    }
    else
#endif
      {
	dip_DLC_energy+=coulomb.Dprefactor*2.*M_PI/volume*(mz*mz);
	fprintf(stderr,"You are not using the P3M method, therefore dp3m.params.epsilon unknown, I assume metallic borders \n");   
      }  
     
    return dip_DLC_energy;
  } else {
    return 0.0;
  }
   
}    
/* ***************************************************************** */

/* -------------------------------------------------------------------------------
   Subroutine to compute the cut-off (NCUT) necessary in the DLC dipolar part
   to get a certain accuracy (acc). We assume particles to have all them a same 
   value of the dipolar momentum modulus (mu_max). mu_max is taken as the largest value of
   mu inside the sytem. If we assum the gap has a width gap_size (within which there is no particles)

   Lz=h+gap_size

   BE CAREFUL:  (1) We assum the short distance for the slab to be in the Z direction
   (2) You must also tune the other 3D method to the same accuracy, otherwise
   it has no sense to have a good accurated result for DLC-dipolar.
 
   ---------------------------------------------------------------------------------- */     

int mdlc_tune(double error)
{
  double de,n,gc,lz,lx,a,fa1,fa2,fa0,h;
  int     kc,limitkc=200,flag;

  MDLC_TRACE(fprintf(stderr, "%d: mdlc_tune().\n", this_node));
 
  n=(double) n_part;
  lz=box_l[2];
  
  a=box_l[0]*box_l[1];
  mpi_bcast_max_mu();   /* we take the maximum dipole in the system, to be sure that the errors in the other case
			   will be equal or less than for this one */
 
  h=dlc_params.h;
  if (h < 0) return ES_ERROR;

  if(h > lz) {
    fprintf(stderr,"tune DLC dipolar: Slab is larger than the box size !!! \n");
    exit(1);
  }
 
  if(fabs(box_l[0]-box_l[1])>0.001) {
    fprintf(stderr,"tune DLC dipolar: box size in x direction is different from y direction !!! \n");
    fprintf(stderr,"The tuning formula requires both to be equal. \n");
    exit(1);
  }

  lx=box_l[0];

  flag=0;
  for(kc=1;kc<limitkc;kc++){
    gc=kc*2.0*PI/lx;
    fa0=sqrt(9.0*exp(+2.*gc*h)*g1_DLC_dip(gc,lz-h)+22.0*g1_DLC_dip(gc,lz)+9.0*exp(-2.0*gc*h)*g1_DLC_dip(gc,lz+h) );
    fa1=0.5*sqrt(PI/(2.0*a))*fa0;
    fa2=g2_DLC_dip(gc,lz);
    de=n*(mu_max*mu_max)/(4.0*(exp(gc*lz)-1.0)) *(fa1+fa2);
    if(de<error) {flag=1;break;}
  }
 
  if(flag==0) {
    fprintf(stderr,"tune DLC dipolar: Sorry, unable to find a proper cut-off for such system and accuracy.\n");
    fprintf(stderr,"Try modifiying the variable limitkc in the c-code: dlc_correction.cpp  ... \n");
    return ES_ERROR;
  }
 
  dlc_params.far_cut=kc;
 
  MDLC_TRACE(fprintf(stderr, "%d: done mdlc_tune().\n", this_node));
 
  return ES_OK;
 
}		

     
//======================================================================================================================
//======================================================================================================================


int mdlc_sanity_checks()
{
#ifdef PARTIAL_PERIODIC  
  char *errtxt;
  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{006 mdlc requires periodicity 1 1 1} ");
    return 1;
  }
#endif  

  // It will be desirable to have a  checking function that check that the slab geometry is such that 
  // the short direction is along the z component.

  return 0;
}
/* ***************************************************************** */

int mdlc_set_params(double maxPWerror, double gap_size, double far_cut)
{
  MDLC_TRACE(fprintf(stderr, "%d: mdlc_set_params().\n", this_node));
  
  dlc_params.maxPWerror = maxPWerror;
  dlc_params.gap_size = gap_size;
  dlc_params.h = box_l[2] - gap_size;
  
  switch (coulomb.Dmethod) {
#ifdef DP3M
  case  DIPOLAR_MDLC_P3M:
  case  DIPOLAR_P3M:
    coulomb.Dmethod =DIPOLAR_MDLC_P3M;
    break;
#endif  
  case  DIPOLAR_MDLC_DS:
  case  DIPOLAR_DS: 
    coulomb.Dmethod =DIPOLAR_MDLC_DS; 
    break;
  default:
    return ES_ERROR;
  }

  dlc_params.far_cut = far_cut;
  if (far_cut != -1) {
    dlc_params.far_calculated = 0;
  }
  else {
    dlc_params.far_calculated = 1;
    if (mdlc_tune(dlc_params.maxPWerror) == ES_ERROR) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{009 mdlc tuning failed, gap size too small} ");
    }
  }
  mpi_bcast_coulomb_params();

  return ES_OK;
}
/* ***************************************************************** */

#endif








     
