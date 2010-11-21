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
/** \file p3m.c  P3M algorithm for long range coulomb interaction.
 *
 *  For more information about the ewald algorithm,
 *  see \ref ewald.h "ewald.h"
*/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "integrate.h"

#include "global.h"
#include "debug.h"
#include "integrate.h"
#include "particle_data.h"
#include "utils.h"
#include "communication.h"
#include "ewald.h"
#include "thermostat.h"
#include "cells.h"
#include "mmm-common.h"

#ifdef ELECTROSTATICS

/************************************************
 * DEFINES
 ************************************************/
   


/************************************************
 * data types
 ************************************************/


/************************************************
 * variables
 ************************************************/

ewald_struct ewald = { 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0 };

/** number of charged particles (only on master node). */
int ewald_sum_qpart=0;
/** Sum of square of charges (only on master node). */
double ewald_sum_q2 = 0.0;
/** square of sum of charges (only on master node). */
double ewald_square_sum_q = 0.0;
/** ewald cache arrays */
/*@{*/
static double *kvec = NULL;
static int *kxfield = NULL;
static int *kyfield = NULL;
static int *kzfield = NULL;
static int total_kvectors = 0;
/*@}*/
/** \name Inverse box dimensions and derived constants */
/*@{*/
static double ux, ux2, uy, uy2, uz;
/*@}*/
/** number of local particles, equals the size of \ref elc::partblk. */
static int n_localpart = 0;

/** structure for storing of sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/** sin/cos caching */ 
/*@{*/
static SCCache **scx = NULL;
static int    n_scx;  
static SCCache **scy = NULL;
static int    n_scy;  
static SCCache **scz = NULL;
static int    n_scz;  
static int    scxoff;
static int    scyoff;
static int    sczoff;
/*@}*/

/** \name ewald sum buffers */
/*@{*/
  double* sums=NULL;
  double* sumc=NULL;
  double* totsums=NULL;
  double* totsumc=NULL;
/*@}*/

/** \name Private Functions */
/************************************************************/
/*@{*/
/** Calculates the kvectors once at the beginning*/
int EWALD_prepare_kfield();
/*@}*/



int EWALD_prepare_kfield() {
  int kymin,kzmin,ksq,kx,ky,kz,totk;
  double rkx,rky,rkz,rksq;
 
/**----------------------------------------------------
       loop over k-vectors. 
       kx ranges over 0 to kmax only.
       ky ranges over 0 to kmax when kx=0 and over
          -kmax to kmax otherwise.
       kz ranges over 1 to kmax when kx=ky=0 and
           over -kmax to kmax otherwise.
--------------------------------------------------------*/

/* measure the length of the array kvec[] */
  totk=0;
  kymin=0;
  kzmin=1;

  //  EWALD_TRACE(fprintf(stderr,"%d: EWALD_prepare_k_field\n",this_node));

  for(kx = 0; kx <= ewald.kmax; kx++) {
    for(ky = kymin; ky <= ewald.kmax; ky++) {
      for(kz = kzmin; kz <= ewald.kmax; kz++) {
        ksq=kx*kx+ky*ky+kz*kz;
        if((ksq < ewald.kmaxsq) && (ksq !=0)) {
          totk++;
        }
      }           
      kzmin = -ewald.kmax;
    }
    kymin=-ewald.kmax;
  }

  kvec  = realloc(kvec,totk*sizeof(double));
  kxfield = realloc(kxfield,totk*sizeof(int));
  kyfield = realloc(kyfield,totk*sizeof(int));
  kzfield = realloc(kzfield,totk*sizeof(int));

 /* calculate the kvec[] */
  totk=0;
  kymin=0;
  kzmin=1;
  for(kx = 0; kx <= ewald.kmax; kx++) {
    rkx = kx/box_l[0];
    for(ky = kymin; ky <= ewald.kmax; ky++) {
      rky = ky/box_l[1];
      for(kz = kzmin; kz <= ewald.kmax; kz++) {
         rkz = kz/box_l[2];
         ksq=kx*kx+ky*ky+kz*kz;
         if((ksq < ewald.kmaxsq) && (ksq !=0)) {
           rksq = rkx*rkx+rky*rky+rkz*rkz;
           kvec[totk] = exp(-rksq*PI*PI/(ewald.alpha*ewald.alpha))/(rksq*box_l[0]*box_l[1]*box_l[2]*2.0*PI);
EWALD_TRACE(fprintf(stderr,"%d: EWALD_prepare_k_field kvec %5i = %18.12g\n",this_node,totk,kvec[totk]));
	   kxfield[totk] = kx;
	   kyfield[totk] = ky;
	   kzfield[totk] = kz;
           totk++;
         }
      }           
      kzmin = -ewald.kmax;
    }
    kymin=-ewald.kmax;
  }

  return totk;
 
}

/************************************************************/

int tclprint_to_result_EWALD(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  double b=ewald.kmax;

  Tcl_PrintDouble(interp, ewald.r_cut, buffer);
  Tcl_AppendResult(interp, "ewald ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, ewald.alpha, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, b, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

int ewald_set_params(double r_cut, double alpha, int kmax)
{
  if(r_cut < 0)
    return -1;

  ewald.r_cut    = r_cut;
  ewald.r_cut_iL = r_cut*box_l_i[0];

  if (alpha > 0) {
    ewald.alpha   = alpha;
    ewald.alpha_L = alpha*box_l[0];
  }
  else
    if (alpha != -1.0)
      return -4;

  if (kmax > 0) {
    ewald.kmax = kmax;
    ewald.kmaxsq = kmax*kmax;
  }
  else return -5;

  mpi_bcast_coulomb_params();

  return 0;
}

/* TODO: This function is not used anywhere. To be removed?  */
int ewald_set_eps(double eps)
{
  ewald.epsilon = eps;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}


int tclcommand_inter_coulomb_parse_ewald(Tcl_Interp * interp, int argc, char ** argv)
{
  double r_cut, alpha;
  int i, kmax;

  coulomb.method = COULOMB_EWALD;
    
#ifdef PARTIAL_PERIODIC
  if(PERIODIC(0) == 0 ||
     PERIODIC(1) == 0 ||
     PERIODIC(2) == 0)
    {
      Tcl_AppendResult(interp, "Need periodicity (1,1,1) with Coulomb EWALD",
		       (char *) NULL);
      return TCL_ERROR;  
    }
#endif

  if (argc < 2) {
    Tcl_AppendResult(interp, "expected: inter coulomb <bjerrum> ewald <r_cut> <alpha> <kmax>",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(! ARG0_IS_D(r_cut))
    return TCL_ERROR;  

  if(argc != 3) {
    Tcl_AppendResult(interp, "wrong # arguments: inter coulomb <bjerrum> ewald <r_cut> <alpha> <kmax>",
		     (char *) NULL);
    return TCL_ERROR;  
  }

  if(! ARG_IS_D(1, alpha))
    return TCL_ERROR;

  if(! ARG_IS_I(2, kmax))
    return TCL_ERROR;

  if ((i = ewald_set_params(r_cut, alpha, kmax)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "r_cut must be positive", (char *) NULL);
      break;
    case -4:
      Tcl_AppendResult(interp, "alpha must be positive", (char *) NULL);
      break;
    case -5:
      Tcl_AppendResult(interp, "kmax must be greater than zero", (char *) NULL);
    default:;
      Tcl_AppendResult(interp, "unspecified error", (char *) NULL);
    }

    return TCL_ERROR;

  }

  return TCL_OK;
}

int EWALD_sanity_checks_boxl() {
  int ret = 0;

  return ret;
}

void EWALD_scaleby_box_l() {
  ewald.r_cut = ewald.r_cut_iL* box_l[0];
  ewald.alpha = ewald.alpha_L * box_l_i[0];
  EWALD_sanity_checks_boxl();
}

int EWALD_sanity_checks()
{
  char *errtxt;
  int ret = 0;

  if (!PERIODIC(0) || !PERIODIC(1) || !PERIODIC(2)) {
    errtxt = runtime_error(128);
    sprintf(errtxt, "{EWALD requires periodicity 1 1 1} ");
    ret = 1;
  }
  
  if( (box_l[0] != box_l[1]) || (box_l[1] != box_l[2]) ) {
    errtxt = runtime_error(128);
    sprintf(errtxt,"{EWALD at present requires a cubic box} ");
    ret = 1;
  }

  if (EWALD_sanity_checks_boxl()) ret = 1;

  return ret;
}

void EWALD_count_charged_particles()
{  
  Cell *cell;
  Particle *part;
  int i,c,np;
  double node_sums[3], tot_sums[3];

  for(i=0;i<3;i++) { node_sums[i]=0.0; tot_sums[i]=0.0;}

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    part = cell->part;
    np   = cell->n;
    for(i=0;i<np;i++) {
      if( part[i].p.q != 0.0 ) {
	node_sums[0] += 1.0;
	node_sums[1] += SQR(part[i].p.q);
	node_sums[2] += part[i].p.q;
      }
    }
  }
  
  MPI_Reduce(node_sums, tot_sums, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  ewald_sum_qpart    = (int)(tot_sums[0]+0.1);
  ewald_sum_q2       = tot_sums[1];
  ewald_square_sum_q = tot_sums[2];
}

void EWALD_init()
{
  ux  = 1/box_l[0];
  ux2 = ux*ux;
  uy  = 1/box_l[1];
  uy2 = uy*uy;  
  uz  = 1/box_l[2];

  if(coulomb.bjerrum == 0.0) {       
    ewald.r_cut    = 0.0;
    ewald.r_cut_iL = 0.0;
    if(this_node==0) 
      EWALD_TRACE(fprintf(stderr,"0: EWALD_init: Bjerrum length is zero.\n");
		fprintf(stderr,"   Electrostatics switched off!\n"));
  }
  else {  
    EWALD_TRACE(fprintf(stderr,"%d: EWALD_init: \n",this_node));

    if (EWALD_sanity_checks()) return;

    EWALD_TRACE(fprintf(stderr,"%d: EWALD_init: preparing kfield\n",this_node));

    total_kvectors=EWALD_prepare_kfield();

    EWALD_TRACE(fprintf(stderr,"%d: EWALD_total_kvectors=%d\n",this_node,total_kvectors));

    EWALD_TRACE(fprintf(stderr,"%d: EWALD initialized\n",this_node));

    EWALD_count_charged_particles();

  }

}

void EWALD_on_resort_particles()
{ 
  int i;

  n_localpart = cells_get_n_particles();

  EWALD_TRACE(fprintf(stderr,"%d: EWALD_on_resort_particles, n_localpart=%d\n",this_node,n_localpart));

  n_scx = (int)(ceil(ewald.kmax) + 1);
  n_scy = (int)(2*ceil(ewald.kmax) + 1);
  n_scz = (int)(2*ceil(ewald.kmax) + 1);
  scxoff = 0;
  scyoff = (int)(ceil(ewald.kmax));
  sczoff = (int)(ceil(ewald.kmax));
  scx=realloc(scx,total_kvectors*sizeof(SCCache*));
  scy=realloc(scy,total_kvectors*sizeof(SCCache*));
  scz=realloc(scz,total_kvectors*sizeof(SCCache*));
  for (i=0; i<total_kvectors; i++) {
    scx[i]=realloc(scx[i],n_localpart*sizeof(SCCache));
    scy[i]=realloc(scy[i],n_localpart*sizeof(SCCache));
    scz[i]=realloc(scz[i],n_localpart*sizeof(SCCache));
  }
  sums= realloc(sums,total_kvectors*sizeof(double));
  sumc= realloc(sumc,total_kvectors*sizeof(double));
  totsums=realloc(totsums,total_kvectors*sizeof(double));
  totsumc=realloc(totsumc,total_kvectors*sizeof(double));
}


double EWALD_calc_kspace_forces(int force_flag, int energy_flag)
{
  Cell *cell;
  Particle *p;
  int i,j,c,np,k,kx,ky,kz;
  int y,z,d_rs=0; /* TODO: d_rs was missing and code couldn't compile. Now it has the arbitrary value of 0, fix it. */
  double rclx, rcly, rclz, sps, spc;
  double fact,tf,tfc,tfs,db_fsum=0 ; /* TODO: db_fsum was missing and code couldn't compile. Now it has the arbitrary value of 0, fix it */ 
  /* Prefactor for force */
  double q,force_prefac;
  /* k space energy */
  double k_space_energy=0.0, node_k_space_energy=0.0;

  EWALD_TRACE(fprintf(stderr,"%d: EWALD_calc_kspace_forces, force flag=%d, energy flag=%d\n",this_node,force_flag,energy_flag));

  rclx=C_2PI/box_l[0];
  rcly=C_2PI/box_l[1];
  rclz=C_2PI/box_l[2];
  
  y=scyoff;
  z=sczoff;

  force_prefac = coulomb.prefactor;

  /* === K Space Calculations === */

  /* === Calculation of k space sums that are common for energy and forces  === */
  if(energy_flag || force_flag) {
    j=0;
    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;
      for(i=0; i<np; i++) {
        scx[  0][j].c= 1.0;
        scx[  0][j].s= 0.0;
        scy[y+0][j].c= 1.0;
        scy[y+0][j].s= 0.0;
        scz[z+0][j].c= 1.0;
        scz[z+0][j].s= 0.0;
        scx[  1][j].c= cos(rclx*p[i].r.p[0]);
        scx[  1][j].s= sin(rclx*p[i].r.p[0]);
        scy[y+1][j].c= cos(rcly*p[i].r.p[1]); 
        scy[y+1][j].s= sin(rcly*p[i].r.p[1]); 
        scz[z+1][j].c= cos(rclz*p[i].r.p[2]); 
        scz[z+1][j].s= sin(rclz*p[i].r.p[2]);
        scy[y-1][j].c= scy[y+1][j].c;
        scy[y-1][j].s=-scy[y+1][j].s;
        scz[z-1][j].c= scz[z+1][j].c;
        scz[z-1][j].s=-scz[z+1][j].s;
        j++;
      }
    }
    

    for (k=2; k<=ewald.kmax; k++) {
      j=0;
      for (c = 0; c < local_cells.n; c++) {
	cell = local_cells.cell[c];
	p  = cell->part;
	np = cell->n;
	for(i=0; i<np; i++) {
          scx[  k][j].c= scx[  k-1][j].c*scx[  1][j].c - scx[  k-1][j].s*scx[  1][j].s;
          scx[  k][j].s= scx[  k-1][j].s*scx[  1][j].c + scx[  k-1][j].c*scx[  1][j].s;
          scy[y+k][j].c= scy[y+k-1][j].c*scy[y+1][j].c - scy[y+k-1][j].s*scy[y+1][j].s;
          scy[y+k][j].s= scy[y+k-1][j].s*scy[y+1][j].c + scy[y+k-1][j].c*scy[y+1][j].s;
          scz[z+k][j].c= scz[z+k-1][j].c*scz[z+1][j].c - scz[z+k-1][j].s*scz[z+1][j].s;
          scz[z+k][j].s= scz[z+k-1][j].s*scz[z+1][j].c + scz[z+k-1][j].c*scz[z+1][j].s;
          scy[y-k][j].c=  scy[y+k][j].c;
          scy[y-k][j].s= -scy[y+k][j].s;
          scz[z-k][j].c=  scz[z+k][j].c;
          scz[z-k][j].s= -scz[z+k][j].s;
          j++;
        }
      }
    }


    for (k=0; k<total_kvectors; k++) {
      kx=(  kxfield[k]);
      ky=(y+kyfield[k]);
      kz=(z+kzfield[k]);
      sums[k]=0.0;
      sumc[k]=0.0;
      j=0;
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i=0; i<np; i++) {
          if (p[i].p.q != 0.0) {
            spc=   scx[kx][j].c*scy[ky][j].c*scz[kz][j].c - scx[kx][j].s*scy[ky][j].s*scz[kz][j].c 
                 - scx[kx][j].c*scy[ky][j].s*scz[kz][j].s - scx[kx][j].s*scy[ky][j].c*scz[kz][j].s;
            sps= - scx[kx][j].s*scy[ky][j].s*scz[kz][j].s + scx[kx][j].c*scy[ky][j].c*scz[kz][j].s
                 + scx[kx][j].c*scy[ky][j].s*scz[kz][j].c + scx[kx][j].s*scy[ky][j].c*scz[kz][j].c;

            sums[k] += p[i].p.q*sps;
            sumc[k] += p[i].p.q*spc;
          }
          j++;
        } 
      }
    }
    MPI_Allreduce(sums,totsums,total_kvectors,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(sumc,totsumc,total_kvectors,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  }


  /* === K Space Energy Calculation  === */
  if(energy_flag) {
    for (k=0; k<total_kvectors; k++) {
      kx=kxfield[k];
      ky=kyfield[k];
      kz=kzfield[k];
      if (kx==0) 
        fact=1.0;
      else
        fact=2.0;
      node_k_space_energy += fact * kvec[k] * (totsums[k]*totsums[k] + totsumc[k]*totsumc[k]);
    }
    EWALD_TRACE(fprintf(stderr,"%d: EWALD: node_k_space_energy=%g\n",this_node,node_k_space_energy));
    node_k_space_energy *= coulomb.prefactor;

    MPI_Reduce(&node_k_space_energy, &k_space_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    EWALD_TRACE(fprintf(stderr,"%d: EWALD: 1 k_space_energy=%g\n",this_node,k_space_energy));

    /* self energy correction */
    k_space_energy -= coulomb.prefactor*(ewald_sum_q2*ewald.alpha_L*box_l_i[0] * wupii);
    EWALD_TRACE(fprintf(stderr,"%d: EWALD: 2 k_space_energy=%g\n",this_node,k_space_energy));
    /* net charge correction */
    /*    k_space_energy -= coulomb.prefactor*(ewald_square_sum_q*PI / (2.0*box_l[0]*SQR(ewald.alpha_L))); */
    EWALD_TRACE(fprintf(stderr,"%d: EWALD: 3 k_space_energy=%g\n",this_node,k_space_energy));
  }
    

  /* === K Space Force Calculation  === */
  if(force_flag) {
    for (k=0; k<total_kvectors; k++) {
      kx=kxfield[k];
      ky=y+kyfield[k];
      kz=z+kzfield[k];
      if (kx==0) 
        fact=1.0;
      else
        fact=2.0;
      tfc=fact*kvec[k]*totsumc[k];
      tfs=fact*kvec[k]*totsums[k];
      j=0;
      for (c = 0; c < local_cells.n; c++) {
	cell = local_cells.cell[c];
	p  = cell->part;
	np = cell->n;
	for(i=0; i<np; i++) { 
	  if( (q=p[i].p.q) != 0.0 ) {
            spc=   scx[kx][j].c*scy[ky][j].c*scz[kz][j].c - scx[kx][j].s*scy[ky][j].s*scz[kz][j].c 
                 - scx[kx][j].s*scy[ky][j].c*scz[kz][j].s - scx[kx][j].c*scy[ky][j].s*scz[kz][j].s;
            sps=   scx[kx][j].s*scy[ky][j].s*scz[kz][j].s - scx[kx][j].c*scy[ky][j].c*scz[kz][j].s
                 - scx[kx][j].c*scy[ky][j].s*scz[kz][j].c - scx[kx][j].s*scy[ky][j].c*scz[kz][j].c;
            tf= spc*tfs+sps*tfc*q;
            p[i].f.f[0] += force_prefac*2.0*C_2PI*tf*kx;
            p[i].f.f[1] += force_prefac*2.0*C_2PI*tf*ky;
            p[i].f.f[2] += force_prefac*2.0*C_2PI*tf*kz;
	  }
          j++;

	    ONEPART_TRACE(if(p[i].p.identity==check_id) fprintf(stderr,"%d: OPT: EWALD  f = (%.3e,%.3e,%.3e) in dir %d add %.5f\n",this_node,p[i].f.f[0],p[i].f.f[1],p[i].f.f[2],d_rs,-db_fsum));

	}
      }
    }
  }

/* currently, only metallic boundary conditions are allowed
  if (ewald.epsilon != EWALD_EPSILON_METALLIC)
    k_space_energy -= calc_dipole_term(force_flag, energy_flag);
*/

  return k_space_energy;
}

/* TODO: This function is not used anywhere. To be removed? */
void   EWALD_exit()
{ 
  /* free memory */

  int i;

  for (i=0; i<total_kvectors; i++) {
    free(scx[i]);
    free(scy[i]);
    free(scz[i]);
  }
  free(scx);
  free(scy);
  free(scz);
  free(sums);
  free(sumc);
  free(totsums);
  free(totsumc);
}

/************************************************************/

/************************************************
 * Debug functions printing ewald structures 
 ************************************************/
/*TODO: this function is not used anywhere. To be removed? */
void ewald_print_ewald_struct(ewald_struct ps) {
  fprintf(stderr,"%d: ewald_struct: \n",this_node);
  fprintf(stderr,"   alpha_L=%f, r_cut_iL=%f \n",
	  ps.alpha_L,ps.r_cut_iL);
}

#endif
