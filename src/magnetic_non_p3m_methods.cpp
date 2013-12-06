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
/** \file magnetic_non_p3m_methods.cpp  All 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas 
 *   Assumes minimum image convention for those axis in which the 
 *   system is periodic as defined by setmd.
 *
 *   MDDS => Calculates dipole-dipole interaction of a perioidic system
 *   by explicitly summing the dipole-dipole interaction over several copies of the system
 *   Uses spherical summation order
 *
 */

#include "domain_decomposition.hpp"
#include "magnetic_non_p3m_methods.hpp"

#ifdef DIPOLES

// Calculates dipolar energy and/or force between two particles
double calc_dipole_dipole_ia(Particle* p1, Particle *p2, int force_flag)
{
  double u,r,pe1,pe2,pe3,pe4,r3,r5,r2,r7,a,b,cc,d,ab;
#ifdef ROTATION
  double bx,by,bz,ax,ay,az; 
#endif
  double  ffx,ffy,ffz;
  double dr[3];
 
	
  // Distance between particles
  get_mi_vector(dr,p1->r.p,p2->r.p);

  // Powers of distance
  r2=dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
  r=sqrt(r2);
  r3=r2*r;
  r5=r3*r2;
  r7=r5*r2;
 
  // Dot products
  pe1=p1->r.dip[0]*p2->r.dip[0]+p1->r.dip[1]*p2->r.dip[1]+p1->r.dip[2]*p2->r.dip[2];
  pe2=p1->r.dip[0]*dr[0]+p1->r.dip[1]*dr[1]+p1->r.dip[2]*dr[2];
  pe3=p2->r.dip[0]*dr[0]+p2->r.dip[1]*dr[1]+p2->r.dip[2]*dr[2];
  pe4=3.0/r5;

  // Energy, if requested
  u= coulomb.Dprefactor* ( pe1/r3 - pe4*pe2*pe3);

  // Force, if requested
  if(force_flag) { 
    a=pe4*pe1;
    b=-15.0*pe2*pe3/r7;
    ab =a+b;
    cc=pe4*pe3;
    d=pe4*pe2;
    
    //  Result
    ffx=ab*dr[0]+cc*p1->r.dip[0]+d*p2->r.dip[0];
    ffy=ab*dr[1]+cc*p1->r.dip[1]+d*p2->r.dip[1];
    ffz=ab*dr[2]+cc*p1->r.dip[2]+d*p2->r.dip[2];
    // Add the force to the particles
    p1->f.f[0] +=coulomb.Dprefactor*ffx;
    p1->f.f[1] +=coulomb.Dprefactor*ffy;
    p1->f.f[2] +=coulomb.Dprefactor*ffz;
    p2->f.f[0] -=coulomb.Dprefactor*ffx;
    p2->f.f[1] -=coulomb.Dprefactor*ffy;
    p2->f.f[2] -=coulomb.Dprefactor*ffz;

    // Torques
#ifdef ROTATION
    ax=p1->r.dip[1]*p2->r.dip[2]-p2->r.dip[1]*p1->r.dip[2];
    ay=p2->r.dip[0]*p1->r.dip[2]-p1->r.dip[0]*p2->r.dip[2];
    az=p1->r.dip[0]*p2->r.dip[1]-p2->r.dip[0]*p1->r.dip[1];
    
    bx=p1->r.dip[1]*dr[2]-dr[1]*p1->r.dip[2];
    by=dr[0]*p1->r.dip[2]-p1->r.dip[0]*dr[2];
    bz=p1->r.dip[0]*dr[1]-dr[0]*p1->r.dip[1];
    
    p1->f.torque[0]+=coulomb.Dprefactor*(-ax/r3+bx*cc);
    p1->f.torque[1]+=coulomb.Dprefactor *(-ay/r3+by*cc);
    p1->f.torque[2]+=coulomb.Dprefactor*(-az/r3+bz*cc);
    
    // 2nd particle     
    bx=p2->r.dip[1]*dr[2]-dr[1]*p2->r.dip[2];
    by=dr[0]*p2->r.dip[2]-p2->r.dip[0]*dr[2];
    bz=p2->r.dip[0]*dr[1]-dr[0]*p2->r.dip[1];
	     
    p2->f.torque[0]+=coulomb.Dprefactor* (ax/r3+bx*d);
    p2->f.torque[1]+=coulomb.Dprefactor*(ay/r3+by*d);
    p2->f.torque[2]+=coulomb.Dprefactor* (az/r3+bz*d);
#endif
  }    
	
  // Return energy
  return u;
}

/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/

double dawaanr_calculations(int force_flag, int energy_flag)
{
  double u; 
  int i,j,c,cc;
  
  if(n_nodes!=1) {fprintf(stderr,"error:  DAWAANR is just for one cpu .... \n"); exit(1);}
  if(!(force_flag) && !(energy_flag) ) {fprintf(stderr," I don't know why you call dawaanr_caclulations with all flags zero \n"); return 0;}
  
  // Variable to sum up the energy
  u=0;

  // Iterate over all cells
  for (c = 0; c < local_cells.n; c++) {
    // Iterate over all particles in this cell
    for(i=0;i<local_cells.cell[c]->n;i++) {
      // If the particle has no dipole moment, ignore it
      if( local_cells.cell[c]->part[i].p.dipm < 1.e-11 ) 
        continue;
      
      // Consider interaction of this particle with the others in the same cell
      for (j=i+1;j<local_cells.cell[c]->n; j++)	{
        // If the particle has no dipole moment, ignore it
        if( local_cells.cell[c]->part[j].p.dipm < 1.e-11 ) 
          continue;
        // Calculate energy and/or force between the particles
	u+=calc_dipole_dipole_ia(&local_cells.cell[c]->part[i],&local_cells.cell[c]->part[j],force_flag);
      }

      // Calculate the ia between this particles and the particles in the 
      // other cells:
      // Iterate over all remaining cells
      for (cc = c+1; cc < local_cells.n; cc++) {
	// Iterate over the particles in this cell
	for (j=0;j<local_cells.cell[cc]->n;j++) {
	  // If it doesn't have dipole moment, ignore
	  if( local_cells.cell[cc]->part[j].p.dipm < 1.e-11 ) 
	    continue;
        
	  // Calculate energy and/or force between the particles
	  u+=calc_dipole_dipole_ia(&local_cells.cell[c]->part[i],&local_cells.cell[cc]->part[j],force_flag);
	}
      }
    }
  }
  
  // Return energy
  return u;
}


/************************************************************/

/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */
/*=================== */

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS               
   =============================================================================
*/

int  Ncut_off_magnetic_dipolar_direct_sum=0;

/************************************************************/


int  magnetic_dipolar_direct_sum_sanity_checks()
{
  /* left for the future , at this moment nothing to do */
  
  return 0;
}

/************************************************************/


double  magnetic_dipolar_direct_sum_calculations(int force_flag, int energy_flag) {
  Cell *cell;
  Particle *part;
  int i,c,np;
  double *x=NULL,  *y=NULL, *z=NULL;
  double *mx=NULL,  *my=NULL, *mz=NULL;
  double *fx=NULL,  *fy=NULL, *fz=NULL;
#ifdef ROTATION
  double *tx=NULL,  *ty=NULL, *tz=NULL;
#endif
  int dip_particles,dip_particles2;
  double ppos[3];
  int img[3];
  double u;

  
  if(n_nodes!=1) {fprintf(stderr,"error: magnetic Direct Sum is just for one cpu .... \n"); exit(1);}
  if(!(force_flag) && !(energy_flag) ) {fprintf(stderr," I don't know why you call dawaanr_caclulations with all flags zero \n"); return 0;}

  x = (double *) malloc(sizeof(double)*n_part);
  y = (double *) malloc(sizeof(double)*n_part);
  z = (double *) malloc(sizeof(double)*n_part);

  mx = (double *) malloc(sizeof(double)*n_part);
  my = (double *) malloc(sizeof(double)*n_part);
  mz = (double *) malloc(sizeof(double)*n_part);
 
  if(force_flag) {
    fx = (double *) malloc(sizeof(double)*n_part);
    fy = (double *) malloc(sizeof(double)*n_part);
    fz = (double *) malloc(sizeof(double)*n_part);
 
#ifdef ROTATION   
    tx = (double *) malloc(sizeof(double)*n_part);
    ty = (double *) malloc(sizeof(double)*n_part);
    tz = (double *) malloc(sizeof(double)*n_part);
#endif  
  }
 
  dip_particles=0;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if( part[i].p.dipm > 1.e-11 ) {
       
	mx[dip_particles]=part[i].r.dip[0];
	my[dip_particles]=part[i].r.dip[1];
	mz[dip_particles]=part[i].r.dip[2];

	/* here we wish the coordinates to be folded into the primary box */
                  
	ppos[0]=part[i].r.p[0];	  
	ppos[1]=part[i].r.p[1];	  
	ppos[2]=part[i].r.p[2];	  
	img[0]=part[i].l.i[0];	  
	img[1]=part[i].l.i[1];	  
        img[2]=part[i].l.i[2];	  		  
        fold_position(ppos, img);
	 
	x[dip_particles]=ppos[0];
	y[dip_particles]=ppos[1];
	z[dip_particles]=ppos[2];
	 

	if(force_flag) {
	  fx[dip_particles]=0;
	  fy[dip_particles]=0;
	  fz[dip_particles]=0;
	 
#ifdef ROTATION
	  tx[dip_particles]=0;
	  ty[dip_particles]=0;
	  tz[dip_particles]=0;
#endif
	} 
	 
	dip_particles++;
	  	 
      }
    }
  }
 
  /*now we do the calculations */
   
  { /* beginning of the area of calculation */
    int nx,ny,nz,i,j;
    double r,rnx,rny,rnz,pe1,pe2,pe3,r3,r5,r2,r7;
    double a,b,c,d;
#ifdef ROTATION
    double ax,ay,az,bx,by,bz;
#endif
    double rx,ry,rz;
    double rnx2,rny2;
    int NCUT[3],NCUT2;
	 
    for(i=0;i<3;i++){
      NCUT[i]=Ncut_off_magnetic_dipolar_direct_sum;
#ifdef PARTIAL_PERIODIC
      if(PERIODIC(i) == 0)  {NCUT[i]=0;}  
#endif            
    }
    NCUT2=Ncut_off_magnetic_dipolar_direct_sum*Ncut_off_magnetic_dipolar_direct_sum;
	     
	   
    u=0;
     
    fprintf(stderr,"Magnetic Direct sum takes long time. Done of %d: \n",dip_particles);

    for(i=0;i<dip_particles;i++){
      fprintf(stderr,"%d\r",i);
      for(j=0;j<dip_particles;j++){
	pe1=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
	rx=x[i]-x[j];
	ry=y[i]-y[j];
	rz=z[i]-z[j];
           
	for(nx=-NCUT[0];nx<=NCUT[0];nx++){
	  rnx=rx+nx*box_l[0]; 
	  rnx2=rnx*rnx;
	  for(ny=-NCUT[1];ny<=NCUT[1];ny++){
	    rny=ry+ny*box_l[1];
	    rny2=rny*rny;
	    for(nz=-NCUT[2];nz<=NCUT[2];nz++){
	      if( !(i==j && nx==0 && ny==0 && nz==0) ) {
		if(nx*nx+ny*ny +nz*nz<=  NCUT2){
		  rnz=rz+nz*box_l[2]; 
		  r2=rnx2+rny2+rnz*rnz;
		  r=sqrt(r2);
		  r3=r2*r;
		  r5=r3*r2;
		  r7=r5*r2;
			    
   
		  pe2=mx[i]*rnx+my[i]*rny+mz[i]*rnz;
		  pe3=mx[j]*rnx+my[j]*rny+mz[j]*rnz;
    
		  /*fprintf(stderr,"--------------------------------\n");
		    fprintf(stderr,"ij: %d %d\n",i,j);
		    fprintf(stderr,"xyz[i]: %lf %lf %lf\n",x[i],y[i],z[i]);
		    fprintf(stderr,"xyz[j]: %lf %lf %lf\n",x[j],y[j],z[j]);
		    fprintf(stderr,"mu xyz[i]: %lf %lf %lf\n",mx[i],my[i],mz[i]);
		    fprintf(stderr,"mu xyz[j]: %lf %lf %lf\n",mx[j],my[j],mz[j]);
		    fprintf(stderr,"rnxyz:  %lf %lf %lf\n",rnx,rny,rnz);
		    fprintf(stderr,"--------------------------------\n");*/
 
    
		  //Energy ............................   
    
		  u+= pe1/r3 - 3.0*pe2*pe3/r5;
	     
		  if(force_flag) {
		    //force ............................
		    a=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
		    a=3.0*a/r5;
		    b=-15.0*pe2*pe3/r7;
		    c=3.0*pe3/r5;
		    d=3.0*pe2/r5;
	     
		    fx[i]+=(a+b)*rnx+c*mx[i]+d*mx[j];
		    fy[i]+=(a+b)*rny+c*my[i]+d*my[j];
		    fz[i]+=(a+b)*rnz+c*mz[i]+d*mz[j];
	     
#ifdef ROTATION
		    //torque ............................
		    c=3.0/r5*pe3;
		    ax=my[i]*mz[j]-my[j]*mz[i];
		    ay=mx[j]*mz[i]-mx[i]*mz[j];
		    az=mx[i]*my[j]-mx[j]*my[i];
	     
		    bx=my[i]*rnz-rny*mz[i];
		    by=rnx*mz[i]-mx[i]*rnz;
		    bz=mx[i]*rny-rnx*my[i];
	     
		    tx[i]+=-ax/r3+bx*c;
		    ty[i]+=-ay/r3+by*c;
		    tz[i]+=-az/r3+bz*c;
#endif	  
		  } /* of force_flag  */
	     
		} }/* of nx*nx+ny*ny +nz*nz< NCUT*NCUT   and   !(i==j && nx==0 && ny==0 && nz==0) */
	    }/* of  for nz */
          }/* of  for ny  */
        }/* of  for nx  */
      }}   /* of  j and i  */ 


    fprintf(stderr,"done \n");
  }/* end of the area of calculation */
    
    
    
  /* set the forces, and torques of the particles within Espresso */
  if(force_flag) {
   
    dip_particles2=0;
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      part = cell->part;
      np   = cell->n;
      for(i=0;i<np;i++) {
	if( part[i].p.dipm  > 1.e-11 ) {
	 
	  part[i].f.f[0]+=coulomb.Dprefactor*fx[dip_particles2];
	  part[i].f.f[1]+=coulomb.Dprefactor*fy[dip_particles2];
	  part[i].f.f[2]+=coulomb.Dprefactor*fz[dip_particles2];
	
#ifdef ROTATION 
	  part[i].f.torque[0]+=coulomb.Dprefactor*tx[dip_particles2];
	  part[i].f.torque[1]+=coulomb.Dprefactor*ty[dip_particles2];
	  part[i].f.torque[2]+=coulomb.Dprefactor*tz[dip_particles2];
#endif
	  dip_particles2++;
	  	 
	}
      }
    }
   
    /* small checking */
    if(dip_particles != dip_particles2) { fprintf(stderr,"magnetic direct sum calculations: error mismatch of particles \n"); exit(1);}
  } /*of if force_flag */
  
  /* free memory used */

  free(x);
  free(y);
  free(z);
  free(mx);
  free(my);
  free(mz);
 
  if(force_flag) {
    free(fx);
    free(fy);
    free(fz);
#ifdef ROTATION
    free(tx);
    free(ty);
    free(tz);
#endif
  }
 
  return 0.5*u;
} 
 
int dawaanr_set_params()
{
  if (n_nodes > 1) {
    return ES_ERROR;
  }
  if (coulomb.Dmethod != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA ) {
    coulomb.Dmethod = DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA;
  } 
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();

  return ES_OK;
}

int mdds_set_params(int n_cut)
{
  if (n_nodes > 1) {
    return ES_ERROR;  
  }
  
  Ncut_off_magnetic_dipolar_direct_sum = n_cut;
  
  if (Ncut_off_magnetic_dipolar_direct_sum == 0) {
    fprintf(stderr,"Careful:  the number of extra replicas to take into account during the direct sum calculation is zero \n");
  }
  
  if (coulomb.Dmethod != DIPOLAR_DS  && coulomb.Dmethod != DIPOLAR_MDLC_DS) {
    coulomb.Dmethod = DIPOLAR_DS;
  }  
  
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
  return ES_OK;
}

#endif
