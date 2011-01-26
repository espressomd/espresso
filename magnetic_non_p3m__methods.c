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
/** \file magnetic_non_p3m__methods.c  All 3d non P3M methods to deal with the magnetic dipoles
 *   
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas exist (i.e. it is not replicated) and we DO NOT 
 *    assume minimum image convention, i.e.,  we work with the unfolded coordinates 
 *
 *   Specific application for which it was intendeed:  handling of a single magnetic filament, where all particles are bonded and hence, it 
 *   it has sense to work with unfolded coordinate, to avoid minimum image convention and replicating the chain when computing the
 *   interactions
 *
 * 
 *  DS => Direct sum , compute the things via direct sum, 
 *
 *  For more information about the DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *  see \ref magnetic_non_p3m__methods.h "magnetic_non_p3m__methods.h"
 */

//#include <mpi.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include "utils.h"
//#include "integrate.h"
//#include "global.h"
//#include "grid.h"
#include "domain_decomposition.h"
//#include "particle_data.h"
//#include "communication.h"


#include "magnetic_non_p3m__methods.h"


#ifdef MAGNETOSTATICS


/* =============================================================================
                  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA                
   =============================================================================
*/
#ifdef DAWAANR

int  tclprint_to_result_DAWAANR(Tcl_Interp *interp){
  Tcl_AppendResult(interp, " dawaanr ", (char *) NULL);
  return TCL_OK;
}

/************************************************************/

int tclcommand_inter_magnetic_parse_dawaanr(Tcl_Interp * interp, int argc, char ** argv)
{
    
  if (coulomb.Dmethod != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA ) {
    coulomb.Dmethod = DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA;
  }  
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with  DAWAANR",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (n_nodes > 1) {
    Tcl_AppendResult(interp, "sorry: DAWAANR only works with 1 cpu", (char *) NULL);
    return TCL_ERROR;  
  }


  coulomb.Dprefactor = (temperature > 0) ? temperature*coulomb.Dbjerrum : coulomb.Dbjerrum;
  return TCL_OK;
}


/************************************************************/


int DAWAANR_sanity_checks()
{
#ifdef PARTIAL_PERIODIC
  char *errtxt;
#endif  
  
  return 0;
}

/************************************************************/

            
double dawaanr_calculations(int force_flag, int energy_flag) { 
  Cell *cell;
  Particle *part;
  int i,c,np;
  double *x=NULL,  *y=NULL, *z=NULL;
  double *mx=NULL,  *my=NULL, *mz=NULL;
  double *fx=NULL,  *fy=NULL, *fz=NULL;
#ifdef ROTATION
  double *tx=NULL,  *ty=NULL, *tz=NULL;
#endif
  int dip_particles,dip_particles2,j;
  double u,r,pe1,pe2,pe3,r3,r5,r2,r7,rx,ry,rz,a,b,cc,d;
#ifdef ROTATION
  double bx,by,bz,ax,ay,az; 
#endif
  double  ffx,ffy,ffz;
 
  if(n_nodes!=1) {fprintf(stderr,"error:  DAWAANR is just for one cpu .... \n"); exit(1);}
  if(!(force_flag) && !(energy_flag) ) {fprintf(stderr," I don't know why you call dawaanr_caclulations with all flags zero \n"); return 0;}

  x = (double *) malloc(sizeof(double)*n_total_particles);
  y = (double *) malloc(sizeof(double)*n_total_particles);
  z = (double *) malloc(sizeof(double)*n_total_particles);
    
  mx = (double *) malloc(sizeof(double)*n_total_particles);
  my = (double *) malloc(sizeof(double)*n_total_particles);
  mz = (double *) malloc(sizeof(double)*n_total_particles);

  if(force_flag) {
 
    fx = (double *) malloc(sizeof(double)*n_total_particles);
    fy = (double *) malloc(sizeof(double)*n_total_particles);
    fz = (double *) malloc(sizeof(double)*n_total_particles);
    
#ifdef ROTATION
    tx = (double *) malloc(sizeof(double)*n_total_particles);
    ty = (double *) malloc(sizeof(double)*n_total_particles);
    tz = (double *) malloc(sizeof(double)*n_total_particles);
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
	 
	/* REMEMBER WE DON'T WANT CALCULATION WITH PERIODIC IMAGES , and we don't have PBC, therefore
	   we must use unfolded positions. Here, in versions prior to 23rd April 2010, there was an error
	   because in fact, we were using only the first part of the rhs do the next 3 expressions, so we where using
	   the almost folded positions.  JJCP MODIFIED on 23rd April 2010 */

	x[dip_particles]=part[i].r.p[0] +part[i].l.i[0]*box_l[0];
	y[dip_particles]=part[i].r.p[1] +part[i].l.i[1]*box_l[1];
	z[dip_particles]=part[i].r.p[2] +part[i].l.i[2]*box_l[2];
	 
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
   
  u=0;
  for(i=0;i<dip_particles-1;i++) {
    for(j=i+1;j<dip_particles;j++) {

      rx=x[i]-x[j];
      ry=y[i]-y[j];
      rz=z[i]-z[j];
 #ifdef PARTIAL_PERIODIC
      // Apply minimum image convention.
      // Due to the layout of the data, we can't use get_mi_vecotr from grid.h
      if (PERIODIC(0))
       rx -= dround(rx*box_l_i[0])*box_l[0];
      if (PERIODIC(1))
       ry -= dround(ry*box_l_i[1])*box_l[1];
      if (PERIODIC(2))
       rz -= dround(rz*box_l_i[2])*box_l[2];
#endif 

      pe1=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
      r2=rx*rx+ry*ry+rz*rz;
	    
      r=sqrt(r2);
      r3=r2*r;
      r5=r3*r2;
      r7=r5*r2;
    
      pe2=mx[i]*rx+my[i]*ry+mz[i]*rz;
      pe3=mx[j]*rx+my[j]*ry+mz[j]*rz;
    
      //Energy ............................   
    
      u+= pe1/r3 - 3.0*pe2*pe3/r5;
	     
      if(force_flag) { 
	//force ............................
	a=mx[i]*mx[j]+my[i]*my[j]+mz[i]*mz[j];
	a=3.0*a/r5;
	b=-15.0*pe2*pe3/r7;
	cc=3.0/r5*pe3;
	d=3.0/r5*pe2;
	     
	// forces  acting on partice i
	ffx=(a+b)*rx+cc*mx[i]+d*mx[j];
	ffy=(a+b)*ry+cc*my[i]+d*my[j];
	ffz=(a+b)*rz+cc*mz[i]+d*mz[j];
	     
	fx[i]+=ffx;
	fy[i]+=ffy;
	fz[i]+=ffz;
	     
	// forces acting on particle j
	fx[j]-=ffx;
	fy[j]-=ffy;
	fz[j]-=ffz;
	     
#ifdef ROTATION
	//torque  acting on particle i ............................

	ax=my[i]*mz[j]-my[j]*mz[i];
	ay=mx[j]*mz[i]-mx[i]*mz[j];
	az=mx[i]*my[j]-mx[j]*my[i];
	     
	bx=my[i]*rz-ry*mz[i];
	by=rx*mz[i]-mx[i]*rz;
	bz=mx[i]*ry-rx*my[i];
	     
	tx[i]+=-ax/r3+bx*cc;
	ty[i]+=-ay/r3+by*cc;
	tz[i]+=-az/r3+bz*cc;
	     
	//torque  acting on particle j ............................
	     
	bx=my[j]*rz-ry*mz[j];
	by=rx*mz[j]-mx[j]*rz;
	bz=mx[j]*ry-rx*my[j];
	     
	tx[j]+=  ax/r3+bx*d;
	ty[j]+= ay/r3+by*d;
	tz[j]+= az/r3+bz*d;
#endif
      } /*of if force_flag */  
    }}

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
    if(dip_particles != dip_particles2) { fprintf(stderr,"dawaanr calculations: error mismatch of particles \n"); exit(1);}
  
  } /* of if force_flag */
  
  
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
  return u;
}
#endif   /*of  ifdef DAWAANR  */

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
#ifdef MAGNETIC_DIPOLAR_DIRECT_SUM

int  Ncut_off_magnetic_dipolar_direct_sum=0;

int tclprint_to_result_Magnetic_dipolar_direct_sum_(Tcl_Interp *interp){
  char buffer[TCL_DOUBLE_SPACE];

  Tcl_AppendResult(interp, " mdds ", buffer, (char *) NULL);
  Tcl_PrintDouble(interp,Ncut_off_magnetic_dipolar_direct_sum , buffer);
  Tcl_AppendResult(interp, " ", buffer, (char *) NULL);

  return TCL_OK;
}

/************************************************************/

int tclcommand_inter_magnetic_parse_mdds(Tcl_Interp * interp, int argc, char ** argv)
{
  int  n_cut=-1;
   
  if (coulomb.Dmethod != DIPOLAR_DS  && coulomb.Dmethod !=DIPOLAR_MDLC_DS ) {
    coulomb.Dmethod = DIPOLAR_DS;
  }  
    

  if (n_nodes > 1) {
    Tcl_AppendResult(interp, "sorry: magnetic dipolar direct sum  only works with 1 cpu", (char *) NULL);
    return TCL_ERROR;  
  }
   
  while(argc > 0) {
    if (ARG0_IS_S("n_cut")) {
      if (! (argc > 1 && ARG1_IS_I(n_cut) && n_cut >= 0)) {
	Tcl_AppendResult(interp, "n_cut expects an nonnegative integer",
			 (char *) NULL);
	return TCL_ERROR;
      } else {
	Ncut_off_magnetic_dipolar_direct_sum=n_cut;
      }    
    } else { /* unknown parameter*/
      Tcl_AppendResult(interp, "unknown parameter/s for the magnetic dipolar direct sum, the only one accepted is:  n_cut  positive_integer", (char *) NULL);
      return TCL_ERROR;
    }
    
    if( Ncut_off_magnetic_dipolar_direct_sum==0) {
      fprintf(stderr,"Careful:  the number of extra replicas to take into account during the direct sum calculation is zero \n");
    }
    
    argc -= 2;
    argv += 2;
  }
   
  coulomb.Dprefactor = (temperature > 0) ? temperature*coulomb.Dbjerrum : coulomb.Dbjerrum; 
  return TCL_OK;
}


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

  x = (double *) malloc(sizeof(double)*n_total_particles);
  y = (double *) malloc(sizeof(double)*n_total_particles);
  z = (double *) malloc(sizeof(double)*n_total_particles);

  mx = (double *) malloc(sizeof(double)*n_total_particles);
  my = (double *) malloc(sizeof(double)*n_total_particles);
  mz = (double *) malloc(sizeof(double)*n_total_particles);
 
  if(force_flag) {
    fx = (double *) malloc(sizeof(double)*n_total_particles);
    fy = (double *) malloc(sizeof(double)*n_total_particles);
    fz = (double *) malloc(sizeof(double)*n_total_particles);
 
#ifdef ROTATION   
    tx = (double *) malloc(sizeof(double)*n_total_particles);
    ty = (double *) malloc(sizeof(double)*n_total_particles);
    tz = (double *) malloc(sizeof(double)*n_total_particles);
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
     

    for(i=0;i<dip_particles;i++){
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
 

#endif /*of ifdef MAGNETIC_DIPOLAR_DIRECT_SUM */




#endif
