/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file mmm2d.cpp  MMM2D algorithm for long range coulomb interaction.
 *
 *  For more information about MMM2D, see \ref mmm2d.hpp "mmm2d.h".
 */

#include <cmath>
#include <mpi.h>
#include "utils.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "cells.hpp"
#include "mmm2d.hpp"
#include "mmm-common.hpp"
#include "specfunc.hpp"
#include "integrate.hpp"
#include "layered.hpp"

#ifdef ELECTROSTATICS
char const *mmm2d_errors[] = {
   "ok",
   "Layer height too large for MMM2D near formula, increase n_layers",
   "box_l[1]/box_l[0] too large for MMM2D near formula, please exchange x and y",
   "Could find not reasonable Bessel cutoff. Please decrease n_layers or the error bound",
   "Could find not reasonable Polygamma cutoff. Consider exchanging x and y",
   "Far cutoff too large, decrease the error bound",
   "Layer height too small for MMM2D far formula, decrease n_layers or skin",
   "IC requires layered cellsystem with more than 3 layers",
};

/** if you define this, the Besselfunctions are calculated up
    to machine precision, otherwise 10^-14, which should be
    definitely enough for daily life. */
#undef BESSEL_MACHINE_PREC

// #define CHECKPOINTS
#if 0
#define LOG_FORCES(x) x
#else
#define LOG_FORCES(x)
#endif

#ifndef BESSEL_MACHINE_PREC
#define K0 LPK0
#define K1 LPK1
#endif

/****************************************
 * LOCAL DEFINES
 ****************************************/

/** Largest reasonable cutoff for far formula. A double cannot overflow
    with this value. */
#define MAXIMAL_FAR_CUT 100

/** Largest reasonable cutoff for Bessel function. The Bessel functions
    are quite slow, so do not make too large. */
#define MAXIMAL_B_CUT 50

/** Largest reasonable order of polygamma series. These are pretty fast,
    so use more of them. Also, the real cutoff is determined at run time,
    so normally we are faster */
#define MAXIMAL_POLYGAMMA 100

/** internal relative precision of far formula. This controls how many
    p,q vectors are done at once. This has nothing to do with the effective
    precision, but rather controls how different values can be we add up without
    loosing the smallest values. In principle one could choose smaller values, but
    that would not make things faster */
#define FARRELPREC 1e-6

/** number of steps in the complex cutoff table */
#define COMPLEX_STEP 16
/** map numbers from 0 to 1/2 onto the complex cutoff table
    (with security margin) */
#define COMPLEX_FAC (COMPLEX_STEP/(.5 + 0.01))

/****************************************
 * LOCAL VARIABLES
 ****************************************/

/** up to that error the sums in the NF are evaluated */
static double part_error;

/** cutoffs for the bessel sum */
static IntList besselCutoff = {NULL, 0, 0};

/** cutoffs for the complex sum */
static int  complexCutoff[COMPLEX_STEP + 1];
/** bernoulli numbers divided by n */
static DoubleList  bon = {NULL, 0, 0};

/** inverse box dimensions */
/*@{*/
static double ux, ux2, uy, uy2, uz;
/*@}*/

/** maximal z for near formula, minimal z for far formula.
    Is identical in the theory, but with the verlet tricks
    this is no longer true, the skin has to be added/subtracted */
/*@{*/
static double max_near, min_far;
/*@}*/

///
static double self_energy;

MMM2D_struct mmm2d_params = { 1e100, 10, 1, 0, 0, 1, 1, 1 };

/** return codes for \ref MMM2D_tune_near and \ref MMM2D_tune_far */
/*@{*/
/** cell too large */
#define ERROR_LARGE 1
/** box too large */
#define ERROR_BOXL 2
/** no reasonable bessel cutoff found */
#define ERROR_BESSEL 3
/** no reasonable polygamma cutoff found */
#define ERROR_POLY 4
/** no reasonable cutoff for the far formula found */
#define ERROR_FARC 5
/** cell too small */
#define ERROR_SMALL 6
/** IC layer requirement */
#define ERROR_ICL 7
/*@}*/

/****************************************
 * LOCAL ARRAYS
 ****************************************/

/** \name Product decomposition data organization
    For the cell blocks
    it is assumed that the lower blocks part is in the lower half.
    This has to have positive sign, so that has to be first. */
/*@{*/

#define POQESP 0
#define POQECP 1
#define POQESM 2
#define POQECM 3

#define PQESSP 0
#define PQESCP 1
#define PQECSP 2
#define PQECCP 3
#define PQESSM 4
#define PQESCM 5
#define PQECSM 6
#define PQECCM 7

#define QQEQQP 0
#define QQEQQM 1

#define ABEQQP 0
#define ABEQZP 1
#define ABEQQM 2
#define ABEQZM 3

/*@}*/

/** number of local particles */
static int n_localpart = 0;

/** temporary buffers for product decomposition */
static double *partblk = NULL;
/** for all local cells including ghosts */
static double *lclcblk = NULL;
/** collected data from the cells above the top neighbor
    of a cell rsp. below the bottom neighbor
    (P=below, M=above, as the signs in the exp). */
static double *gblcblk = NULL;

/** contribution from the image charges */
static double lclimge[8]; 

typedef struct {
  double s, c;
} SCCache;

/** sin/cos caching */ 
static SCCache *scxcache = NULL;
static int    n_scxcache;  
/** sin/cos caching */ 
static SCCache *scycache = NULL;
static int    n_scycache;  


/** \name Local functions for the near formula */
/************************************************************/
/*@{*/

/** complex evaluation */
static void prepareBernoulliNumbers(int nmax);

/** cutoff error setup. Returns error code */
static int MMM2D_tune_near(double error);

/** energy of all local particles with their copies */
void MMM2D_self_energy();

/*@}*/

/** \name Local functions for the far formula */
/************************************************************/
/*@{*/

/** sin/cos storage */
static void prepare_scx_cache();
static void prepare_scy_cache();
/** clear the image contributions if there is no dielectric contrast and no image charges */
static void clear_image_contributions(int size);
/** gather the informations for the far away image charges */
static void gather_image_contributions(int size);
/** spread the top/bottom sums */
static void distribute(int size, double fac);
/** 2 pi |z| code */
static void setup_z_force();
static void setup_z_energy();
static void   add_z_force();
static double     z_energy();
/** p=0 per frequency code */
static void setup_P(int p, double omega, double fac);
static void   add_P_force();
static double     P_energy(double omega);
/** q=0 per frequency code */
static void setup_Q(int q, double omega, double fac);
static void   add_Q_force();
static double     Q_energy(double omega);
/** p,q <> 0 per frequency code */
static void setup_PQ(int p, int q, double omega, double fac);
static void   add_PQ_force(int p, int q, double omega);
static double     PQ_energy(double omega);

/** cutoff error setup. Returns error code */
static int MMM2D_tune_far(double error);

/*@}*/

///
void MMM2D_setup_constants()
{
  ux  = 1/box_l[0];
  ux2 = ux*ux;
  uy  = 1/box_l[1];
  uy2 = uy*uy;  
  uz  = 1/box_l[2];

  switch (cell_structure.type) {
  case CELL_STRUCTURE_NSQUARE:
    max_near = box_l[2];
    /* not used */
    min_far = 0.0;
    break;
  case CELL_STRUCTURE_LAYERED:
    max_near = 2*layer_h + skin;
    min_far  =   layer_h - skin;    
    break;
  default:
    fprintf(stderr, "%d: INTERNAL ERROR: MMM2D setup for cell structure it should reject\n", this_node);
    errexit();
  }
}

/****************************************
 * FAR FORMULA
 ****************************************/

static void prepare_scx_cache()
{
  int np, c, i, ic, freq, o;
  double pref, arg;
  Particle *part;
  
  for (freq = 1; freq <= n_scxcache; freq++) {
    pref = C_2PI*ux*freq;
    o = (freq-1)*n_localpart;
    ic = 0;
    for (c = 1; c <= n_layers; c++) {
      np   = cells[c].n;
      part = cells[c].part;
      for (i = 0; i < np; i++) {
	arg = pref*part[i].r.p[0];
	scxcache[o + ic].s = sin(arg);
	scxcache[o + ic].c = cos(arg);
	ic++;
      }
    }
  }
}

static void prepare_scy_cache()
{
  int np, c, i, ic, freq, o;
  double pref, arg;
  Particle *part;
  
  for (freq = 1; freq <= n_scycache; freq++) {
    pref = C_2PI*uy*freq;
    o = (freq-1)*n_localpart;
    ic = 0;
    for (c = 1; c <= n_layers; c++) {
      np   = cells[c].n;
      part = cells[c].part;
      for (i = 0; i < np; i++) {
	arg = pref*part[i].r.p[1];
	scycache[o + ic].s = sin(arg);
	scycache[o + ic].c = cos(arg);
	ic++;
      }
    }
  }
}

/*****************************************************************/
/* data distribution */
/*****************************************************************/

/* vector operations */

/** pdc = 0 */
inline void clear_vec(double *pdc, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc[i] = 0;
}

/** pdc_d = pdc_s */
inline void copy_vec(double *pdc_d, double *pdc_s, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

/** pdc_d = pdc_s1 + pdc_s2 */
inline void add_vec(double *pdc_d, double *pdc_s1, double *pdc_s2, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc_s1 + pdc_s2 */
inline void addscale_vec(double *pdc_d, double scale, double *pdc_s1, double *pdc_s2, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = scale*pdc_s1[i] + pdc_s2[i];
}

/** pdc_d = scale*pdc */
inline void scale_vec(double scale, double *pdc, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc[i] *= scale;
}

/* block indexing - has to fit to the PQ block definitions above.
   size gives the full size of the data block,
   e_size is the size of only the top or bottom half, i.e. half of size.
*/

inline double *block(double *p, int index, int size)
{
  return &p[index*size];
}

inline double *blwentry(double *p, int index, int e_size)
{
  return &p[2*index*e_size];
}

inline double *abventry(double *p, int index, int e_size)
{
  return &p[(2*index + 1)*e_size];
}

/* dealing with the image contributions from far outside the simulation box */

void clear_image_contributions(int e_size)
{
  if (this_node == 0)
    /* the gblcblk contains all contributions from layers deeper than one layer below our system,
       which is precisely what the gblcblk should contain for the lowest layer. */
    clear_vec(blwentry(gblcblk, 0, e_size), e_size);

  if (this_node == n_nodes - 1)
    /* same for the top node */
    clear_vec(abventry(gblcblk, n_layers - 1, e_size), e_size);
}

void gather_image_contributions(int e_size)
{
  double recvbuf[8];

  /* collect the image charge contributions with at least a layer distance */
  MPI_Allreduce(lclimge, recvbuf, 2*e_size, MPI_DOUBLE, MPI_SUM, comm_cart);

  if (this_node == 0)
    /* the gblcblk contains all contributions from layers deeper than one layer below our system,
       which is precisely what the gblcblk should contain for the lowest layer. */
    copy_vec(blwentry(gblcblk, 0, e_size),recvbuf, e_size);

  if (this_node == n_nodes - 1)
    /* same for the top node */
    copy_vec(abventry(gblcblk, n_layers - 1, e_size), recvbuf + e_size, e_size);
}

/* the data transfer routine for the lclcblks itself */
void distribute(int e_size, double fac)
{
  int c, node, inv_node;
  double sendbuf[8];
  double recvbuf[8];
  MPI_Status status;

  /* send/recv to/from other nodes. Also builds up the gblcblk. */
  for (node = 0; node < n_nodes; node++) {
    inv_node = n_nodes - node - 1;
    /* up */
    if (node == this_node) {
      /* calculate sums of cells below */
      for (c = 1; c < n_layers; c++)
	addscale_vec(blwentry(gblcblk, c, e_size), fac, blwentry(gblcblk, c - 1, e_size), blwentry(lclcblk, c - 1, e_size), e_size);

      /* calculate my ghost contribution only if a node above exists */
      if (node + 1 < n_nodes) {
	addscale_vec(sendbuf, fac, blwentry(gblcblk, n_layers - 1, e_size), blwentry(lclcblk, n_layers - 1, e_size), e_size);
	copy_vec(sendbuf + e_size, blwentry(lclcblk, n_layers, e_size), e_size);
	MPI_Send(sendbuf, 2*e_size, MPI_DOUBLE, node + 1, 0, comm_cart);
      }
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*e_size, MPI_DOUBLE, node, 0, comm_cart, &status);
      copy_vec(blwentry(gblcblk, 0, e_size), recvbuf, e_size);
      copy_vec(blwentry(lclcblk, 0, e_size), recvbuf + e_size, e_size);
    }

    /* down */
    if (inv_node == this_node) {
      /* calculate sums of all cells above */
      for (c = n_layers + 1; c > 2; c--)
	addscale_vec(abventry(gblcblk, c - 3, e_size), fac, abventry(gblcblk, c - 2, e_size), abventry(lclcblk, c, e_size), e_size);
      
      /* calculate my ghost contribution only if a node below exists */
      if (inv_node -  1 >= 0) {
	addscale_vec(sendbuf, fac, abventry(gblcblk, 0, e_size), abventry(lclcblk, 2, e_size), e_size);
	copy_vec(sendbuf + e_size, abventry(lclcblk, 1, e_size), e_size);
	MPI_Send(sendbuf, 2*e_size, MPI_DOUBLE, inv_node - 1, 0, comm_cart);
      }
    }
    else if (inv_node - 1 == this_node) {
      MPI_Recv(recvbuf, 2*e_size, MPI_DOUBLE, inv_node, 0, comm_cart, &status);
      copy_vec(abventry(gblcblk, n_layers - 1, e_size), recvbuf, e_size);
      copy_vec(abventry(lclcblk, n_layers + 1, e_size), recvbuf + e_size, e_size);
    }
  }
}

#ifdef CHECKPOINTS
static void checkpoint(char *text, int p, int q, int e_size)
{
  int c, i;
  fprintf(stderr, "%d: %s %d %d\n", this_node, text, p, q);

  fprintf(stderr, "partblk\n");
  for (c = 0; c < n_localpart; c++) {
    fprintf(stderr, "%d", c);    
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(partblk, c, 2*e_size)[i]);
    fprintf(stderr, " m");
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(partblk, c, 2*e_size)[i + e_size]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "lclcblk\n");
  fprintf(stderr, "0");    
  for (i = 0; i < e_size; i++)
    fprintf(stderr, " %10.3g", block(lclcblk, 0, 2*e_size)[i]);
  fprintf(stderr, "\n");
  for (c = 1; c <= n_layers; c++) {
    fprintf(stderr, "%d", c);    
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(lclcblk, c, 2*e_size)[i]);
    fprintf(stderr, " m");
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(lclcblk, c, 2*e_size)[i + e_size]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "%d", n_layers + 1);
  for (i = 0; i < e_size; i++)
    fprintf(stderr, "           ");
  fprintf(stderr, " m");

  for (i = 0; i < e_size; i++)
    fprintf(stderr, " %10.3g", block(lclcblk, n_layers + 1, 2*e_size)[i + e_size]);
  fprintf(stderr, "\n");

  fprintf(stderr, "gblcblk\n");
  for (c = 0; c < n_layers; c++) {
    fprintf(stderr, "%d", c + 1);    
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(gblcblk, c, 2*e_size)[i]);
    fprintf(stderr, " m");
    for (i = 0; i < e_size; i++)
      fprintf(stderr, " %10.3g", block(gblcblk, c, 2*e_size)[i + e_size]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}
#else
#define checkpoint(text,p,q,size)
#endif

/*****************************************************************/
/* 2 pi (sign)(z) */
/*****************************************************************/

static void setup_z_force()
{
  int np, c, i;
  double pref = coulomb.prefactor*C_2PI*ux*uy;
  Particle *part;
  double *lclimgebot=NULL,*lclimgetop=NULL;
  int e_size=1,size = 2;
  double e, e_di_l, e_di_h;

  double fac_imgsum;

  /* in case of metallic boundary conditions on both sides, we get an infinite array,
     which only exists for charge neutral systems. But in this case, we can as well
     not sum up the force array, as the net force per image is 0 */
  if (mmm2d_params.delta_mult != 1.0) {
    fac_imgsum = 1/(1 - mmm2d_params.delta_mult);
  }
  else {
    fac_imgsum = 0;
  }

  if (mmm2d_params.dielectric_contrast_on) 
    clear_vec(lclimge, size); 

  if(this_node==0) {
    
    lclimgebot=blwentry(lclcblk,0,e_size);
    clear_vec(lclimgebot, e_size);
  }
  
  if(this_node==n_nodes-1) {
    
    lclimgetop=abventry(lclcblk,n_layers+1,e_size);
    clear_vec(lclimgetop, e_size);
  }

  /* calculate local cellblks. partblks don't make sense */
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    lclcblk[size*c] = 0;
    for (i = 0; i < np; i++) {
      lclcblk[size*c] += part[i].p.q;
      
      if (mmm2d_params.dielectric_contrast_on) {
	e_di_l = (mmm2d_params.delta_mult*mmm2d_params.delta_mid_bot
		  + mmm2d_params.delta_mult)*fac_imgsum;
	if (c==1 && this_node==0) {
	  e = mmm2d_params.delta_mid_bot;  
	  lclimgebot[QQEQQP] += part[i].p.q*e;
	}
	else
	  e_di_l += mmm2d_params.delta_mid_bot;

	e_di_h = (mmm2d_params.delta_mult*mmm2d_params.delta_mid_top
		  + mmm2d_params.delta_mult)*fac_imgsum;

	if (c==n_layers && this_node==n_nodes-1) {
	  e = mmm2d_params.delta_mid_top;
	  lclimgetop[QQEQQP] += part[i].p.q*e;
	}
	else
	  e_di_h += mmm2d_params.delta_mid_top;

	lclimge[QQEQQP] += part[i].p.q*e_di_l;
       	lclimge[QQEQQM] += part[i].p.q*e_di_h;
      }   
    }
    lclcblk[size*c] *= pref;
    lclcblk[size*c+1] = lclcblk[size*c];   
  }

  if (mmm2d_params.dielectric_contrast_on) {
    scale_vec(pref, lclimge, size);
    if(this_node==0)
      scale_vec(pref, blwentry(lclcblk, 0, e_size), e_size);
    if(this_node==n_nodes-1)
      scale_vec(pref, abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
}

static void add_z_force()
{
  int np, c, i;
  double add;
  Particle *part;
  double *othcblk;
  int size = 2;

  for (c = 1; c <= n_layers; c++) {
    othcblk = block(gblcblk, c - 1, size);
    add = othcblk[QQEQQP] - othcblk[QQEQQM];
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      part[i].f.f[2] += part[i].p.q*add;
      LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
			 this_node, part[i].p.identity, part[i].f.f[0],
			 part[i].f.f[1], part[i].f.f[2]));
    }
  }
}

static void setup_z_energy()
{
  int np, c, i;
  double pref = -coulomb.prefactor*C_2PI*ux*uy;
  Particle *part;
  int e_size = 2, size = 4;

  if (this_node == 0)
    /* the lowest lclcblk does not contain anything, since there are no charges below the simulation
       box, at least for this term. */
    clear_vec(blwentry(lclcblk, 0, e_size), e_size);

  if (this_node == n_nodes - 1)
    /* same for the top node */
    clear_vec(abventry(lclcblk, n_layers + 1, e_size), e_size);

  /* calculate local cellblks. partblks don't make sense */
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    clear_vec(blwentry(lclcblk, c, e_size), e_size);
    for (i = 0; i < np; i++) {
      lclcblk[size*c + ABEQQP] += part[i].p.q;
      lclcblk[size*c + ABEQZP] += part[i].p.q*part[i].r.p[2];
    }
    scale_vec(pref, blwentry(lclcblk, c, e_size), e_size);
    /* just to be able to use the standard distribution. Here below
       and above terms are the same */
    copy_vec(abventry(lclcblk, c, e_size), blwentry(lclcblk, c, e_size), e_size);
  }
}

static double z_energy()
{
  int np, c, i;
  Particle *part;
  double *othcblk;
  int size = 4;
  double eng = 0;
  for (c = 1; c <= n_layers; c++) {
    othcblk = block(gblcblk, c - 1, size);
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      eng += part[i].p.q*(part[i].r.p[2]*othcblk[ABEQQP] - othcblk[ABEQZP] -
			  part[i].r.p[2]*othcblk[ABEQQM] + othcblk[ABEQZM]);
      
    }
  }
  return eng;
}

/*****************************************************************/
/* PoQ exp sum */
/*****************************************************************/
static void setup_P(int p, double omega, double fac)
{
  int np, c, i, ic, o = (p-1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*4*M_PI*ux*uy*fac*fac;
  double h = box_l[2];
  double fac_imgsum = 1/(1 - mmm2d_params.delta_mult*exp(-omega*2*h));
  double fac_delta_mid_bot = mmm2d_params.delta_mid_bot*fac_imgsum; 
  double fac_delta_mid_top = mmm2d_params.delta_mid_top*fac_imgsum;
  double fac_delta         = mmm2d_params.delta_mult*fac_imgsum;
  double layer_top;
  double e, e_di_l, e_di_h;
  double *llclcblk;
  double *lclimgebot = NULL, *lclimgetop = NULL;
  int e_size = 2, size = 4;

  if (mmm2d_params.dielectric_contrast_on)
    clear_vec(lclimge, size); 

  if(this_node==0) {
    /* on the lowest node, clear the lclcblk below, which only contains the images of the lowest layer
       if there is dielectric contrast, otherwise it is empty */
    lclimgebot = block(lclcblk, 0, size);
    clear_vec(blwentry(lclcblk, 0, e_size), e_size);
  }
  if(this_node==n_nodes-1) {
    /* same for the top node */
    lclimgetop = block(lclcblk, n_layers + 1, size);
    clear_vec(abventry(lclcblk, n_layers + 1, e_size), e_size);
  }

  layer_top = my_left[2] + layer_h;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - layer_top));

      partblk[size*ic + POQESM] = part[i].p.q*scxcache[o + ic].s/e;
      partblk[size*ic + POQESP] = part[i].p.q*scxcache[o + ic].s*e;
      partblk[size*ic + POQECM] = part[i].p.q*scxcache[o + ic].c/e;
      partblk[size*ic + POQECP] = part[i].p.q*scxcache[o + ic].c*e;

      /* take images due to different dielectric constants into account */
      if (mmm2d_params.dielectric_contrast_on) {
	if (c==1 && this_node==0) {
	  /* There are image charges at -(2h+z) and -(2h-z) etc. layer_h included due to the shift
	     in z */
	  e_di_l = ( exp(omega*(-part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_bot +
		     exp(omega*( part[i].r.p[2] - 2*h + layer_h))                   )*fac_delta;

	  e = exp(omega*(-part[i].r.p[2]))*mmm2d_params.delta_mid_bot;

	  lclimgebot[POQESP] += part[i].p.q*scxcache[o + ic].s*e;
	  lclimgebot[POQECP] += part[i].p.q*scxcache[o + ic].c*e;
	}
	else
	  /* There are image charges at -(z) and -(2h-z) etc. layer_h included due to the shift in z */
	  e_di_l = ( exp(omega*(-part[i].r.p[2]       + layer_h)) +
		     exp(omega*( part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_top )*fac_delta_mid_bot;    

	if (c==n_layers && this_node==n_nodes-1) {
	  /* There are image charges at (3h-z) and (h+z) from the top layer etc. layer_h included
	     due to the shift in z */
	  e_di_h = ( exp(omega*( part[i].r.p[2] - 3*h + 2*layer_h))*mmm2d_params.delta_mid_top +
		     exp(omega*(-part[i].r.p[2] -   h + 2*layer_h))                 )*fac_delta;
	  
	  /* There are image charges at (h-z) layer_h included due to the shift in z */
	  e = exp(omega*(part[i].r.p[2] - h + layer_h))*mmm2d_params.delta_mid_top;
	  
	  lclimgetop[POQESM]+= part[i].p.q*scxcache[o + ic].s*e;
	  lclimgetop[POQECM]+= part[i].p.q*scxcache[o + ic].c*e;
	}
	else
	  /* There are image charges at (h-z) and (h+z) from the top layer etc. layer_h included due
	     to the shift in z */
	  e_di_h = ( exp(omega*( part[i].r.p[2] - h + 2*layer_h)) +
		     exp(omega*(-part[i].r.p[2] - h + 2*layer_h))*mmm2d_params.delta_mid_bot )*fac_delta_mid_top;

	lclimge[POQESP] += part[i].p.q*scxcache[o + ic].s*e_di_l;
	lclimge[POQECP] += part[i].p.q*scxcache[o + ic].c*e_di_l;
	lclimge[POQESM] += part[i].p.q*scxcache[o + ic].s*e_di_h;
	lclimge[POQECM] += part[i].p.q*scxcache[o + ic].c*e_di_h;
      }

      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, blwentry(lclcblk, c, e_size), e_size);
    scale_vec(pref, abventry(lclcblk, c, e_size), e_size);

    layer_top += layer_h;
  }

  if (mmm2d_params.dielectric_contrast_on) {
    scale_vec(pref, lclimge, size);
    if(this_node==0)
      scale_vec(pref, blwentry(lclcblk, 0, e_size), e_size);
    if(this_node==n_nodes-1)
      scale_vec(pref, abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
}

/* compare setup_P */
static void setup_Q(int q, double omega, double fac)
{
  int np, c, i, ic, o = (q-1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*4*M_PI*ux*uy*fac*fac;
  double h = box_l[2];
  double fac_imgsum = 1/(1 - mmm2d_params.delta_mult*exp(-omega*2*h));
  double fac_delta_mid_bot = mmm2d_params.delta_mid_bot*fac_imgsum; 
  double fac_delta_mid_top = mmm2d_params.delta_mid_top*fac_imgsum;
  double fac_delta         = mmm2d_params.delta_mult*fac_imgsum;
  double layer_top;
  double e, e_di_l, e_di_h;
  double *llclcblk;
  double *lclimgebot=NULL, *lclimgetop=NULL;
  int e_size = 2, size = 4;
 
  if (mmm2d_params.dielectric_contrast_on)
    clear_vec(lclimge, size); 
  
  if(this_node==0) {
    lclimgebot = block(lclcblk, 0, size);
    clear_vec(blwentry(lclcblk, 0, e_size), e_size);
  }
  
  if(this_node==n_nodes-1) {
    lclimgetop = block(lclcblk, n_layers + 1, size);
    clear_vec(abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
  
  layer_top = my_left[2] + layer_h;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - layer_top));

      partblk[size*ic + POQESM] = part[i].p.q*scycache[o + ic].s/e;
      partblk[size*ic + POQESP] = part[i].p.q*scycache[o + ic].s*e;
      partblk[size*ic + POQECM] = part[i].p.q*scycache[o + ic].c/e;
      partblk[size*ic + POQECP] = part[i].p.q*scycache[o + ic].c*e;

      if (mmm2d_params.dielectric_contrast_on) {
	if(c==1 && this_node==0) {
	  e_di_l = ( exp(omega*(-part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_bot +
		     exp(omega*( part[i].r.p[2] - 2*h + layer_h))                 )*fac_delta;
	
	  e = exp(omega*(-part[i].r.p[2]))*mmm2d_params.delta_mid_bot;

	  lclimgebot[POQESP] += part[i].p.q*scycache[o + ic].s*e;
	  lclimgebot[POQECP] += part[i].p.q*scycache[o + ic].c*e;
	}
	else
	  e_di_l = ( exp(omega*(-part[i].r.p[2]       + layer_h)) +
		     exp(omega*( part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_top )*fac_delta_mid_bot;    

	if(c==n_layers && this_node==n_nodes-1) {
	  e_di_h = (exp(omega*( part[i].r.p[2] - 3*h + 2*layer_h))*mmm2d_params.delta_mid_top +
		    exp(omega*(-part[i].r.p[2] -   h + 2*layer_h))                 )*fac_delta;
	  	  
	  e = exp(omega*(part[i].r.p[2] - h + layer_h))*mmm2d_params.delta_mid_top;
	  
	  lclimgetop[POQESM] += part[i].p.q*scycache[o + ic].s*e;
	  lclimgetop[POQECM] += part[i].p.q*scycache[o + ic].c*e;
	}
	else
	  e_di_h = ( exp(omega*( part[i].r.p[2] - h + 2*layer_h)) +
		     exp(omega*(-part[i].r.p[2] - h + 2*layer_h))*mmm2d_params.delta_mid_bot )*fac_delta_mid_top;

	lclimge[POQESP] += part[i].p.q*scycache[o + ic].s*e_di_l;
	lclimge[POQECP] += part[i].p.q*scycache[o + ic].c*e_di_l;
	lclimge[POQESM] += part[i].p.q*scycache[o + ic].s*e_di_h;
	lclimge[POQECM] += part[i].p.q*scycache[o + ic].c*e_di_h;
      }
      
      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, blwentry(lclcblk, c, e_size), e_size);
    scale_vec(pref, abventry(lclcblk, c, e_size), e_size);

    layer_top += layer_h;
  }

  if (mmm2d_params.dielectric_contrast_on) {
    scale_vec(pref, lclimge, size);
    if(this_node==0)
      scale_vec(pref, blwentry(lclcblk, 0, e_size), e_size);
    if(this_node==n_nodes-1)
      scale_vec(pref, abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
}

static void add_P_force()
{
  int np, c, i, ic;
  Particle *part;
  double *othcblk;
  int size = 4;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    othcblk = block(gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
      part[i].f.f[0] +=
	partblk[size*ic + POQESM]*othcblk[POQECP] - partblk[size*ic + POQECM]*othcblk[POQESP] +
	partblk[size*ic + POQESP]*othcblk[POQECM] - partblk[size*ic + POQECP]*othcblk[POQESM];
      part[i].f.f[2] +=
	partblk[size*ic + POQECM]*othcblk[POQECP] + partblk[size*ic + POQESM]*othcblk[POQESP] -
	partblk[size*ic + POQECP]*othcblk[POQECM] - partblk[size*ic + POQESP]*othcblk[POQESM];

      LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
			 this_node, part[i].p.identity, part[i].f.f[0],
			 part[i].f.f[1], part[i].f.f[2]));
      ic++;
    }
  }
}

static double P_energy(double omega)
{
  int np, c, i, ic;
  double *othcblk;
  int size = 4;
  double eng = 0;
  double pref = 1/omega;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    othcblk = block(gblcblk, c - 1, size);
    for (i = 0; i < np; i++) {
      eng += pref*(partblk[size*ic + POQECM]*othcblk[POQECP] + partblk[size*ic + POQESM]*othcblk[POQESP] +
		   partblk[size*ic + POQECP]*othcblk[POQECM] + partblk[size*ic + POQESP]*othcblk[POQESM]);
      ic++;
    }
  }
  return eng;
}

static void add_Q_force()
{
  int np, c, i, ic;
  Particle *part;
  double *othcblk;
  int size = 4;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    othcblk = block(gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
      part[i].f.f[1] +=
	partblk[size*ic + POQESM]*othcblk[POQECP] - partblk[size*ic + POQECM]*othcblk[POQESP] +
	partblk[size*ic + POQESP]*othcblk[POQECM] - partblk[size*ic + POQECP]*othcblk[POQESM];
      part[i].f.f[2] +=
	partblk[size*ic + POQECM]*othcblk[POQECP] + partblk[size*ic + POQESM]*othcblk[POQESP] -
	partblk[size*ic + POQECP]*othcblk[POQECM] - partblk[size*ic + POQESP]*othcblk[POQESM];

      LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
			 this_node, part[i].p.identity, part[i].f.f[0],
			 part[i].f.f[1], part[i].f.f[2]));
      ic++;
    }
  }
}

static double Q_energy(double omega)
{
  int np, c, i, ic;
  double *othcblk;
  int size = 4;
  double eng = 0;
  double pref = 1/omega;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    othcblk = block(gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
      eng += pref*(partblk[size*ic + POQECM]*othcblk[POQECP] + partblk[size*ic + POQESM]*othcblk[POQESP] +
		   partblk[size*ic + POQECP]*othcblk[POQECM] + partblk[size*ic + POQESP]*othcblk[POQESM]);
      ic++;
    }
  }
  return eng;
}

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

/* compare setup_P */
static void setup_PQ(int p, int q, double omega, double fac)
{
  int np, c, i, ic, ox = (p - 1)*n_localpart, oy = (q - 1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*8*M_PI*ux*uy*fac*fac;
  double h = box_l[2];
  double fac_imgsum = 1/(1 - mmm2d_params.delta_mult*exp(-omega*2*h));
  double fac_delta_mid_bot = mmm2d_params.delta_mid_bot*fac_imgsum; 
  double fac_delta_mid_top = mmm2d_params.delta_mid_top*fac_imgsum;
  double fac_delta         = mmm2d_params.delta_mult*fac_imgsum;
  double layer_top;
  double e, e_di_l, e_di_h;
  double *llclcblk;
  double *lclimgebot=NULL, *lclimgetop=NULL;
  int e_size = 4, size = 8;

  if (mmm2d_params.dielectric_contrast_on)
    clear_vec(lclimge, size); 

  if(this_node==0) {
    lclimgebot = block(lclcblk, 0, size);
    clear_vec(blwentry(lclcblk, 0, e_size), e_size);
  }
  
  if(this_node==n_nodes-1) {
    lclimgetop = block(lclcblk, n_layers + 1, size);
    clear_vec(abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
  
  layer_top = my_left[2] + layer_h;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - layer_top));
      
      partblk[size*ic + PQESSM] = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q/e;
      partblk[size*ic + PQESCM] = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q/e;
      partblk[size*ic + PQECSM] = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q/e;
      partblk[size*ic + PQECCM] = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q/e;

      partblk[size*ic + PQESSP] = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e;
      partblk[size*ic + PQESCP] = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e;
      partblk[size*ic + PQECSP] = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e;
      partblk[size*ic + PQECCP] = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e;

      if (mmm2d_params.dielectric_contrast_on) {
	if(c==1 && this_node==0) {	
	  e_di_l = (exp(omega*(-part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_bot +
		    exp(omega*( part[i].r.p[2] - 2*h + layer_h))                 )*fac_delta;

	  e = exp(omega*(-part[i].r.p[2]))*mmm2d_params.delta_mid_bot;
	
	  lclimgebot[PQESSP] += scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e;
	  lclimgebot[PQESCP] += scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e;
	  lclimgebot[PQECSP] += scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e;
	  lclimgebot[PQECCP] += scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e;
	}
	else	
	  e_di_l = ( exp(omega*(-part[i].r.p[2]       + layer_h)) +
		     exp(omega*( part[i].r.p[2] - 2*h + layer_h))*mmm2d_params.delta_mid_top )*fac_delta_mid_bot;     
	
	if(c==n_layers && this_node==n_nodes-1) {
	  e_di_h = ( exp(omega*( part[i].r.p[2] - 3*h + 2*layer_h))*mmm2d_params.delta_mid_top +
		     exp(omega*(-part[i].r.p[2] -   h + 2*layer_h))                 )*fac_delta;
	  
	  e = exp(omega*(part[i].r.p[2]-h+layer_h))*mmm2d_params.delta_mid_top;
	
	  lclimgetop[PQESSM] += scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e;
	  lclimgetop[PQESCM] += scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e;
	  lclimgetop[PQECSM] += scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e;
	  lclimgetop[PQECCM] += scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e;
	}
	else
	  e_di_h = ( exp(omega*( part[i].r.p[2] - h + 2*layer_h)) +
		     exp(omega*(-part[i].r.p[2] - h + 2*layer_h))*mmm2d_params.delta_mid_bot )*fac_delta_mid_top;  
      
        lclimge[PQESSP] += scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e_di_l;
	lclimge[PQESCP] += scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e_di_l;
	lclimge[PQECSP] += scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e_di_l;
	lclimge[PQECCP] += scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e_di_l;

        lclimge[PQESSM] += scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e_di_h;
	lclimge[PQESCM] += scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e_di_h;
	lclimge[PQECSM] += scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e_di_h;
	lclimge[PQECCM] += scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e_di_h;
      }
      
      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, blwentry(lclcblk, c, e_size), e_size);
    scale_vec(pref, abventry(lclcblk, c, e_size), e_size);
    
    layer_top += layer_h;
  }

  if (mmm2d_params.dielectric_contrast_on) {
    scale_vec(pref, lclimge, size);

    if(this_node==0)
      scale_vec(pref, blwentry(lclcblk, 0, e_size), e_size);
    if(this_node==n_nodes-1)
      scale_vec(pref, abventry(lclcblk, n_layers + 1, e_size), e_size);
  }
}

static void add_PQ_force(int p, int q, double omega)
{
  int np, c, i, ic;
  Particle *part;
  double pref_x = C_2PI*ux*p/omega;
  double pref_y = C_2PI*uy*q/omega;
  double *othcblk;
  int size = 8;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    othcblk = block(gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
      part[i].f.f[0] +=
	pref_x*(partblk[size*ic + PQESCM]*othcblk[PQECCP] + partblk[size*ic + PQESSM]*othcblk[PQECSP] -
		partblk[size*ic + PQECCM]*othcblk[PQESCP] - partblk[size*ic + PQECSM]*othcblk[PQESSP] +
		partblk[size*ic + PQESCP]*othcblk[PQECCM] + partblk[size*ic + PQESSP]*othcblk[PQECSM] -
		partblk[size*ic + PQECCP]*othcblk[PQESCM] - partblk[size*ic + PQECSP]*othcblk[PQESSM]);
      part[i].f.f[1] +=
	pref_y*(partblk[size*ic + PQECSM]*othcblk[PQECCP] + partblk[size*ic + PQESSM]*othcblk[PQESCP] -
		partblk[size*ic + PQECCM]*othcblk[PQECSP] - partblk[size*ic + PQESCM]*othcblk[PQESSP] +
		partblk[size*ic + PQECSP]*othcblk[PQECCM] + partblk[size*ic + PQESSP]*othcblk[PQESCM] -
		partblk[size*ic + PQECCP]*othcblk[PQECSM] - partblk[size*ic + PQESCP]*othcblk[PQESSM]);
      part[i].f.f[2] +=
	       (partblk[size*ic + PQECCM]*othcblk[PQECCP] + partblk[size*ic + PQECSM]*othcblk[PQECSP] +
	        partblk[size*ic + PQESCM]*othcblk[PQESCP] + partblk[size*ic + PQESSM]*othcblk[PQESSP] -
	        partblk[size*ic + PQECCP]*othcblk[PQECCM] - partblk[size*ic + PQECSP]*othcblk[PQECSM] -
	        partblk[size*ic + PQESCP]*othcblk[PQESCM] - partblk[size*ic + PQESSP]*othcblk[PQESSM]);

      LOG_FORCES(fprintf(stderr, "%d: part %d force %10.3g %10.3g %10.3g\n",
			 this_node, part[i].p.identity, part[i].f.f[0],
			 part[i].f.f[1], part[i].f.f[2]));
      ic++;
    }
  }
}

static double PQ_energy(double omega)
{
  int size = 8;
  double eng = 0;
  double pref = 1/omega;

  int ic = 0;
  for (int c = 1; c <= n_layers; c++) {
    int np   = cells[c].n;
    double *othcblk = block(gblcblk, c - 1, size);

    for (int i = 0; i < np; i++) {
      eng += pref*(partblk[size*ic + PQECCM]*othcblk[PQECCP] + partblk[size*ic + PQECSM]*othcblk[PQECSP] +
		   partblk[size*ic + PQESCM]*othcblk[PQESCP] + partblk[size*ic + PQESSM]*othcblk[PQESSP] +
		   partblk[size*ic + PQECCP]*othcblk[PQECCM] + partblk[size*ic + PQECSP]*othcblk[PQECSM] +
		   partblk[size*ic + PQESCP]*othcblk[PQESCM] + partblk[size*ic + PQESSP]*othcblk[PQESSM]);
      ic++;
    }
  }
  return eng;
}

/*****************************************************************/
/* main loops */
/*****************************************************************/

static void add_force_contribution(int p, int q)
{
  double omega, fac;
  
  if (q == 0) {
    if (p == 0) {
      setup_z_force();

      if (mmm2d_params.dielectric_contrast_on)
	gather_image_contributions(1);
      else
	clear_image_contributions(1);

      distribute(1, 1.);
      add_z_force();
      checkpoint("************2piz", 0, 0, 1);
    }
    else {
      omega = C_2PI*ux*p;
      fac = exp(-omega*layer_h);
      setup_P(p, omega, fac);
      if (mmm2d_params.dielectric_contrast_on)
	gather_image_contributions(2);
      else
	clear_image_contributions(2);
      distribute(2, fac); 
      add_P_force();
      checkpoint("************distri p", p, 0, 2);
    }
  }
  else if (p == 0) {
    omega = C_2PI*uy*q;
    fac = exp(-omega*layer_h);
    setup_Q(q, omega, fac);
    if (mmm2d_params.dielectric_contrast_on)
      gather_image_contributions(2);
    else
      clear_image_contributions(2);
    distribute(2, fac);
    add_Q_force();
    checkpoint("************distri q", 0, q, 2);
  }
  else {
    omega = C_2PI*sqrt(SQR(ux*p) + SQR(uy*q));
    fac = exp(-omega*layer_h);
    setup_PQ(p, q, omega, fac);
    if (mmm2d_params.dielectric_contrast_on)
      gather_image_contributions(4);
    else
      clear_image_contributions(4);
    distribute(4, fac);
    add_PQ_force(p, q, omega);
    checkpoint("************distri pq", p, q, 4);
  }
}

static double energy_contribution(int p, int q)
{
  double eng;
  double omega, fac;
  
  if (q == 0) {
    if (p == 0) {
      setup_z_energy();
      clear_image_contributions(2);
      distribute(2, 1.);
      eng = z_energy();
      checkpoint("E************2piz", 0, 0, 2);
    }
    else {
      omega = C_2PI*ux*p;
      fac = exp(-omega*layer_h);
      setup_P(p, omega, fac);
      if (mmm2d_params.dielectric_contrast_on)
	gather_image_contributions(2);
      else
	clear_image_contributions(2);
      distribute(2, fac);
      eng = P_energy(omega);
      checkpoint("************distri p", p, 0, 2);
    }
  }
  else if (p == 0) {
    omega = C_2PI*uy*q;
    fac = exp(-omega*layer_h);
    setup_Q(q, omega, fac);
    if (mmm2d_params.dielectric_contrast_on)
      gather_image_contributions(2);
    else
      clear_image_contributions(2);
    distribute(2, fac);
    eng = Q_energy(omega);
    checkpoint("************distri q", 0, q, 2);
  }
  else {
    omega = C_2PI*sqrt(SQR(ux*p) + SQR(uy*q));
    fac = exp(-omega*layer_h);
    setup_PQ(p, q, omega, fac);
    if (mmm2d_params.dielectric_contrast_on)
      gather_image_contributions(4);
    else
      clear_image_contributions(4);
    distribute(4, fac);
    eng = PQ_energy(omega);
    checkpoint("************distri pq", p, q, 4);
  }
  return eng;
}

double MMM2D_add_far(int f, int e)
{
  double eng;
  int p, q;
  double R, dR, q2;
  int *undone;
  
  // It's not really far...
  eng = e ? self_energy : 0;

  if (mmm2d_params.far_cut == 0.0)
    return 0.5*eng;
  
  undone = (int*)malloc((n_scxcache + 1)*sizeof(int));

  prepare_scx_cache();
  prepare_scy_cache();

  /* complicated loop. We work through the p,q vectors in rings
     from outside to inside to avoid problems with cancellation */

  /* up to which q vector we have to work */
  for (p = 0; p <= n_scxcache; p++) {
    if (p == 0)
      q =  n_scycache;
    else {
      q2 = mmm2d_params.far_cut2 - SQR(ux*(p - 1));
      if (q2 > 0)
	q = 1 + box_l[1]*(int)ceil(sqrt(q2));
      else
	q = 1;
      /* just to be on the safe side... */
      if (q > n_scycache) q = n_scycache;
    }
    undone[p] = q;
  }

  dR = -log(FARRELPREC)/C_2PI*uz;
  
  for(R = mmm2d_params.far_cut; R > 0; R -= dR) {
    for (p = n_scxcache; p >= 0; p--) {
      for (q = undone[p]; q >= 0; q--) {
	if (ux2*SQR(p)  + uy2*SQR(q) < SQR(R))
	  break;
	if (f)
	  add_force_contribution(p, q);
	if (e)
	  eng += energy_contribution(p, q);
      }
      undone[p] = q;
    }
  }
  // printf("yyyy\n");
  /* clean up left overs */
  for (p = n_scxcache; p >= 0; p--) {
    q = undone[p];
    // fprintf(stderr, "left over %d\n", q);
    for (; q >= 0; q--) {
      // printf("xxxxx %d %d\n", p, q);
      if (f)
	add_force_contribution(p, q);
      if (e)
	eng += energy_contribution(p, q);
     }
  }
  
  free(undone);
  
  return 0.5*eng;
}

static int MMM2D_tune_far(double error)
{
  double err;
  double min_inv_boxl = dmin(ux, uy);
  mmm2d_params.far_cut = min_inv_boxl;
  do {
    err = exp(-2*M_PI*mmm2d_params.far_cut*min_far)/min_far*
      (C_2PI*mmm2d_params.far_cut + 2*(ux + uy) + 1/min_far);
    // fprintf(stderr, "far tuning: %f -> %f, %f\n", mmm2d_params.far_cut, err, min_far);
    mmm2d_params.far_cut += min_inv_boxl;
  }
  while (err > error && mmm2d_params.far_cut*layer_h < MAXIMAL_FAR_CUT);
  if (mmm2d_params.far_cut*layer_h >= MAXIMAL_FAR_CUT)
    return ERROR_FARC;
  // fprintf(stderr, "far cutoff %g %g %g\n", mmm2d_params.far_cut, err, min_far);
  mmm2d_params.far_cut -= min_inv_boxl;
  mmm2d_params.far_cut2 = SQR(mmm2d_params.far_cut);
  return 0;
}

/****************************************
 * NEAR FORMULA
 ****************************************/

static int MMM2D_tune_near(double error)
{
  int P, n, i;
  double uxrho2m2max, uxrhomax2;
  int p;
  double T, pref, err, exponent;
  double L, sum;

  /* yes, it's y only... */
  if (max_near > box_l[1]/2)
    return ERROR_LARGE;
  if (min_far < 0)
    return ERROR_SMALL;
  if (ux*box_l[1] >= 3/M_SQRT2 )
    return ERROR_BOXL;

  /* error is split into three parts:
     one part for bessel, one for complex
     and one for polygamma cutoff */
  part_error = error/3;

  /* Bessel sum, determine cutoff */
  P = 2;
  exponent = M_PI*ux*box_l[1];
  T  = exp(exponent)/exponent;
  pref = 8*ux*dmax(C_2PI*ux, 1);
  do {
    L = M_PI*ux*(P - 1);
    sum = 0;
    for (p = 1; p <= P; p++)
      sum += p*exp(-exponent*p);
    err = pref*K1(box_l[1]*L)*(T*((L + uy)/M_PI*box_l[0] - 1) + sum);
    P++;
  }
  while (err > part_error && (P - 1) < MAXIMAL_B_CUT);
  P--;
  if (P == MAXIMAL_B_CUT)
    return ERROR_BESSEL;

  // fprintf(stderr, "bessel cutoff %d %g\n", P, err);

  realloc_intlist(&besselCutoff, besselCutoff.n = P);
  for (p = 1; p < P; p++)
    besselCutoff.e[p-1] = (int)floor(((double)P)/(2*p)) + 1;

  /* complex sum, determine cutoffs (dist dependent) */
  T = log(part_error/(16*M_SQRT2)*box_l[0]*box_l[1]);
  // for 0, the sum is exactly zero, so do not calculate anything
  complexCutoff[0] = 0;
  for (i = 1; i <= COMPLEX_STEP; i++)
    complexCutoff[i] = (int)ceil(T/log(i/COMPLEX_FAC));
  prepareBernoulliNumbers(complexCutoff[COMPLEX_STEP]);

  /* polygamma, determine order */
  n = 1;
  uxrhomax2 = SQR(ux*box_l[1])/2;
  uxrho2m2max = 1.0;
  do {
    create_mod_psi_up_to(n+1);

    err = 2*n*fabs(mod_psi_even(n, 0.5))*uxrho2m2max;
    uxrho2m2max *= uxrhomax2;
    n++;
  }
  while (err > 0.1*part_error && n < MAXIMAL_POLYGAMMA);
  if (n == MAXIMAL_POLYGAMMA)
    return ERROR_POLY;
  // fprintf(stderr, "polygamma cutoff %d %g\n", n, err);

  return 0;
}

static void prepareBernoulliNumbers(int bon_order)
{
  int l;
  /* BernoulliB[2 n]/(2 n)!(2 Pi)^(2n) up to order 33 */
  static double bon_table[34] = {
     1.0000000000000000000, 3.2898681336964528729,
    -2.1646464674222763830, 2.0346861239688982794,
    -2.0081547123958886788, 2.0019891502556361707,
    -2.0004921731066160966, 2.0001224962701174097,
    -2.0000305645188173037, 2.0000076345865299997,
    -2.0000019079240677456, 2.0000004769010054555,
    -2.0000001192163781025, 2.0000000298031096567,
    -2.0000000074506680496, 2.0000000018626548648,
    -2.0000000004656623667, 2.0000000001164154418,
    -2.0000000000291038438, 2.0000000000072759591,
    -2.0000000000018189896, 2.0000000000004547474,
    -2.0000000000001136868, 2.0000000000000284217,
    -2.0000000000000071054, 2.0000000000000017764,
    -2.0000000000000004441, 2.0000000000000001110,
    -2.0000000000000000278, 2.0000000000000000069,
    -2.0000000000000000017, 2.0000000000000000004,
    -2.0000000000000000001, 2.0000000000000000000
  };

  if (bon_order < 2)
    bon_order = 2;

  realloc_doublelist(&bon, bon.n = bon_order);

  /* the ux is multiplied in to bessel, complex and psi at once, not here,
     and we use uy*(z + iy), so the uy is also treated below */
  for(l = 1; (l <= bon_order) && (l < 34); l++)
    bon.e[l-1] = 2*uy*bon_table[l];

  for (; l <= bon_order; l++) {
    if (l & 1)
      bon.e[l-1] =  4.0*uy;
    else
      bon.e[l-1] = -4.0*uy;      
  }
}

void add_mmm2d_coulomb_pair_force(double charge_factor,
				  double d[3], double dl2, double dl, double force[3])
{
  double F[3];
  double pref = coulomb.prefactor*charge_factor;
  double z2   = d[2]*d[2];
  double rho2 = d[1]*d[1] + z2;
  int i;

  if (d[2] >box_l[1]/2) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{024 near formula called for too distant particle pair} ");
    return;
  }

  F[0] = F[1] = F[2] = 0;

  /* Bessel sum */
  {
    int p, l;
    double k0, k1;
    double k0Sum, k1ySum, k1Sum;
    double freq;
    double rho_l, ypl;
    double c, s;

    for (p = 1; p < besselCutoff.n; p++) {
      k0Sum  = 0;
      k1ySum = 0;
      k1Sum  = 0;

      freq = C_2PI*ux*p;

      for (l = 1; l < besselCutoff.e[p-1]; l++) {
	ypl   = d[1] + l*box_l[1];
	rho_l = sqrt(ypl*ypl + z2);
#ifdef BESSEL_MACHINE_PREC
	k0 = K0(freq*rho_l);
	k1 = K1(freq*rho_l);
#else
	LPK01(freq*rho_l, &k0, &k1);
#endif
	k1 /= rho_l;
	k0Sum  += k0;
	k1Sum  += k1;
	k1ySum += k1*ypl;

	ypl   = d[1] - l*box_l[1];
	rho_l = sqrt(ypl*ypl + z2);
#ifdef BESSEL_MACHINE_PREC
	k0 = K0(freq*rho_l);
	k1 = K1(freq*rho_l);
#else
	LPK01(freq*rho_l, &k0, &k1);
#endif
	k1 /= rho_l;
	k0Sum  += k0;
	k1Sum  += k1;
	k1ySum += k1*ypl;
      }

      /* the ux is multiplied in to bessel, complex and psi at once, not here */
      c = 4*freq*cos(freq*d[0]);
      s = 4*freq*sin(freq*d[0]);
      F[0] +=      s*k0Sum;
      F[1] +=      c*k1ySum;
      F[2] += d[2]*c*k1Sum;
    }
    // fprintf(stderr, " bessel force %f %f %f\n", F[0], F[1], F[2]);
  }

  /* complex sum */
  {
    double zeta_r, zeta_i;
    double zet2_r, zet2_i;
    double ztn_r,  ztn_i;
    double tmp_r;
    int end, n;

    ztn_r = zeta_r = uy*d[2];
    ztn_i = zeta_i = uy*d[1];
    zet2_r = zeta_r*zeta_r - zeta_i*zeta_i;
    zet2_i = 2*zeta_r*zeta_i;

    end = (int)ceil(COMPLEX_FAC*uy2*rho2);
    if (end > COMPLEX_STEP) {
      end = COMPLEX_STEP;
      fprintf(stderr, "MMM2D: some particles left the assumed slab, precision might be lost\n");
    }
    if (end < 0) {
      char *errtxt = runtime_error(100);
      ERROR_SPRINTF(errtxt, "{MMM2D: distance was negative, coordinates probably out of range } ");
      end = 0;
    }
    end = complexCutoff[end];

    for (n = 0; n < end; n++) {
      F[1] -= bon.e[n]*ztn_i;
      F[2] += bon.e[n]*ztn_r;

      tmp_r = ztn_r*zet2_r - ztn_i*zet2_i;
      ztn_i = ztn_r*zet2_i + ztn_i*zet2_r;
      ztn_r = tmp_r;
    }
    // fprintf(stderr, "complex force %f %f %f %d\n", F[0], F[1], F[2], end);
  }

  /* psi sum */
  {
    int n;
    double uxx = ux*d[0];
    double uxrho2 = ux2*rho2;
    double uxrho_2n, uxrho_2nm2; /* rho^{2n-2} */
    double mpe, mpo;

    /* n = 0 inflicts only Fx and pot */
    /* one ux is multiplied in to bessel, complex and psi at once, not here */
    F[0] += ux*mod_psi_odd(0, uxx);

    uxrho_2nm2 = 1.0;
    for (n = 1;n < n_modPsi; n++) {
      mpe    = mod_psi_even(n, uxx);
      mpo    = mod_psi_odd(n, uxx);
      uxrho_2n = uxrho_2nm2*uxrho2;

      F[0] +=     ux *uxrho_2n  *mpo;
      F[1] += 2*n*ux2*uxrho_2nm2*mpe*d[1];
      F[2] += 2*n*ux2*uxrho_2nm2*mpe*d[2];

      /* y < rho => ux2*uxrho_2nm2*d[1] < ux*uxrho_2n */
      if (fabs(2*n*ux*uxrho_2n*mpe) < part_error)
	break;

      uxrho_2nm2 = uxrho_2n;
    }
    // fprintf(stderr, "    psi force %f %f %f %d\n", F[0], F[1], F[2], n);
  }


  for (i = 0; i < 3; i++)
    F[i] *= ux;

  /* explicitly added potentials r_{-1,0} and r_{1,0} */
  {
    double cx    = d[0] + box_l[0];
    double rinv2 = 1.0/(cx*cx + rho2), rinv = sqrt(rinv2);
    double rinv3 = rinv*rinv2;
    F[0] +=   cx*rinv3;
    F[1] += d[1]*rinv3;
    F[2] += d[2]*rinv3;

    cx   = d[0] - box_l[0];
    rinv2 = 1.0/(cx*cx + rho2); rinv = sqrt(rinv2);
    rinv3 = rinv*rinv2;
    F[0] +=   cx*rinv3;
    F[1] += d[1]*rinv3;
    F[2] += d[2]*rinv3;

    rinv3 = 1/(dl2*dl);
    F[0] += d[0]*rinv3;
    F[1] += d[1]*rinv3;
    F[2] += d[2]*rinv3;

    // fprintf(stderr, "explcit force %f %f %f\n", F[0], F[1], F[2]);
  }

  for (i = 0; i < 3; i++)
    force[i] += pref*F[i];
}

inline double calc_mmm2d_copy_pair_energy(double d[3])
{
  double eng;
  double z2     = d[2]*d[2];
  double rho2   = d[1]*d[1] + z2;

  /* the ux is multiplied in below */
  eng = -2*log(4*M_PI*uy*box_l[0]);

  /* Bessel sum */
  {
    int p, l;
    double k0Sum;
    double freq;
    double rho_l, ypl;
    double c;

    for (p = 1; p < besselCutoff.n; p++) {
      k0Sum  = 0;

      freq = C_2PI*ux*p;

      for (l = 1; l < besselCutoff.e[p-1]; l++) {
	ypl   = d[1] + l*box_l[1];
	rho_l = sqrt(ypl*ypl + z2);
	k0Sum  += K0(freq*rho_l);

	ypl   = d[1] - l*box_l[1];
	rho_l = sqrt(ypl*ypl + z2);
	k0Sum  += K0(freq*rho_l);
      }

      /* the ux is multiplied in to bessel, complex and psi at once, not here */
      c = 4*cos(freq*d[0]);
      eng += c*k0Sum;
    }
    // fprintf(stderr, " bessel energy %f\n", eng);
  }

  /* complex sum */
  {
    double zeta_r, zeta_i;
    double zet2_r, zet2_i;
    double ztn_r,  ztn_i;
    double tmp_r;
    int end, n;

    ztn_r = zeta_r = uy*d[2];
    ztn_i = zeta_i = uy*d[1];

    zet2_r = zeta_r*zeta_r - zeta_i*zeta_i;
    zet2_i = 2*zeta_r*zeta_i;

    ztn_r = zet2_r;
    ztn_i = zet2_i;

    end = (int)ceil(COMPLEX_FAC*uy2*rho2);
    if (end > COMPLEX_STEP) {
	end = COMPLEX_STEP;
	fprintf(stderr, "MMM2D: some particles left the assumed slab, precision might be lost\n");
    }
    end = complexCutoff[end];
    for (n = 1; n <= end; n++) {
      eng -= box_l[1]/(2*n)*bon.e[n-1]*ztn_r;
 
      tmp_r = ztn_r*zet2_r - ztn_i*zet2_i;
      ztn_i = ztn_r*zet2_i + ztn_i*zet2_r;
      ztn_r = tmp_r;
    }
    // fprintf(stderr, "complex energy %f %d\n", eng, end);
  }

  /* psi sum */
  {
    int n;
    double add;
    double uxx = ux*d[0];
    double uxrho2 = ux2*rho2;
    double uxrho_2n;

    /* n = 0 inflicts only Fx and pot */
    /* one ux is multiplied in to bessel, complex and psi at once, not here */
    eng -= mod_psi_even(0, uxx);

    uxrho_2n = uxrho2;
    for (n = 1; n < n_modPsi; n++) {
      add = uxrho_2n*mod_psi_even(n, uxx);
      eng -= add;
      if (fabs(add) < part_error)
	break;
      uxrho_2n *= uxrho2;
    }
    // fprintf(stderr, "    psi energy %f %d\n", eng, n);
  }

  eng *= ux;

  /* explicitly added potentials r_{-1,0} and r_{1,0} */
  {
    double cx   = d[0] + box_l[0];
    double rinv = sqrt(1.0/(cx*cx + rho2));
    eng += rinv;

    cx   = d[0] - box_l[0];
    rinv = sqrt(1.0/(cx*cx + rho2));
    eng += rinv;

    // fprintf(stderr, "explcit energy %f %f %f %f\n", d[0], d[1], d[2], eng);
  }

  return eng;
}

double mmm2d_coulomb_pair_energy(double charge_factor,
				 double dv[3], double d2, double d)
{
  double eng, pref = coulomb.prefactor*charge_factor;
  if (pref != 0.0) {
    eng = calc_mmm2d_copy_pair_energy(dv);
    return pref*(eng + 1/d);
  }
  return 0.0;
}

void MMM2D_self_energy()
{
  int c, np, i;
  Particle *part;
  double dv[3] = {0, 0, 0};
  double seng = coulomb.prefactor*calc_mmm2d_copy_pair_energy(dv);
  self_energy = 0;

  /* this one gives twice the real self energy, as it is used
     in the far formula which counts everything twice and in
     the end divides by two*/

  // fprintf(stderr, "%d: self energy %g\n", this_node, seng);

  for (c = 0; c < local_cells.n; c++) {
    np   = local_cells.cell[c]->n;
    part = local_cells.cell[c]->part;

    for (i = 0; i < np; i++) {
      self_energy += seng*SQR(part[i].p.q);
     }
  }
}

/****************************************
 * COMMON PARTS
 ****************************************/

int MMM2D_set_params(double maxPWerror, double far_cut, double delta_top, double delta_bot)
{
  int err;

  if (cell_structure.type != CELL_STRUCTURE_NSQUARE &&
      cell_structure.type != CELL_STRUCTURE_LAYERED) {
    return ES_ERROR;
  }

  mmm2d_params.maxPWerror = maxPWerror;
    
  if (delta_top != 0.0 || delta_bot != 0.0) {
    mmm2d_params.dielectric_contrast_on = 1;
    mmm2d_params.delta_mid_top = delta_top;
    mmm2d_params.delta_mid_bot = delta_bot;
    mmm2d_params.delta_mult = delta_top*delta_bot;
  }
  else {
    mmm2d_params.dielectric_contrast_on = 0;
    mmm2d_params.delta_mid_top = 0;
    mmm2d_params.delta_mid_bot = 0;
    mmm2d_params.delta_mult = 0;
  }

  MMM2D_setup_constants();

  if ((err = MMM2D_tune_near(maxPWerror)))
    return err;

  /* if we cannot do the far formula, force off */
  if (cell_structure.type == CELL_STRUCTURE_NSQUARE ||
      (cell_structure.type == CELL_STRUCTURE_LAYERED && n_nodes*n_layers < 3)) {
    mmm2d_params.far_cut = 0.0;
    if (mmm2d_params.dielectric_contrast_on) {
      return ERROR_ICL;
    }
  }
  else {
    mmm2d_params.far_cut = far_cut;
    mmm2d_params.far_cut2 = SQR(far_cut);
    if (mmm2d_params.far_cut > 0)
      mmm2d_params.far_calculated = 0;
    else {
      if ((err = MMM2D_tune_far(maxPWerror)))
	return err;
      mmm2d_params.far_calculated = 1;    
    }
  }

  coulomb.method = COULOMB_MMM2D;

  mpi_bcast_coulomb_params();

  return ES_OK;
}

int MMM2D_sanity_checks()
{
  char *errtxt;

  if (!PERIODIC(0) || !PERIODIC(1) || PERIODIC(2)) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{025 MMM2D requires periodicity 1 1 0} ");
    return 1;
  }
  
  if (cell_structure.type != CELL_STRUCTURE_LAYERED &&
      cell_structure.type != CELL_STRUCTURE_NSQUARE) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{026 MMM2D at present requires layered (or n-square) cellsystem} ");
    return 1;
  }
  return 0;
}

void MMM2D_init()
{
  int err;
  char *errtxt;

  if (MMM2D_sanity_checks()) return;

  MMM2D_setup_constants();
  if ((err = MMM2D_tune_near(mmm2d_params.maxPWerror))) {
    errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{027 MMM2D auto-retuning: %s} ", mmm2d_errors[err]);
    coulomb.method = COULOMB_NONE;
    return;
  }
  if (cell_structure.type == CELL_STRUCTURE_NSQUARE ||
      (cell_structure.type == CELL_STRUCTURE_LAYERED && n_nodes*n_layers < 3)) {
    mmm2d_params.far_cut = 0.0;
    if (mmm2d_params.dielectric_contrast_on) {
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{027 MMM2D auto-retuning: IC requires layered cellsystem with > 3 layers} ");
    }
  }
  else {
    if (mmm2d_params.far_calculated) {
      if ((err = MMM2D_tune_far(mmm2d_params.maxPWerror))) {
	errtxt = runtime_error(128);
	ERROR_SPRINTF(errtxt, "{028 MMM2D auto-retuning: %s} ", mmm2d_errors[err]);
	coulomb.method = COULOMB_NONE;
	return;
      }
    }
  }
}

void MMM2D_on_resort_particles()
{
  /* if we need MMM2D far formula, allocate caches */
  if (cell_structure.type == CELL_STRUCTURE_LAYERED) {
    n_localpart = cells_get_n_particles();
    n_scxcache = (int)(ceil(mmm2d_params.far_cut/ux) + 1);
    n_scycache = (int)(ceil(mmm2d_params.far_cut/uy) + 1);
    scxcache = (SCCache*)realloc(scxcache, n_scxcache*n_localpart*sizeof(SCCache));
    scycache = (SCCache*)realloc(scycache, n_scycache*n_localpart*sizeof(SCCache));
    
    partblk   = (double*)realloc(partblk,  n_localpart*8*sizeof(double));
    lclcblk   = (double*)realloc(lclcblk,  n_cells*8*sizeof(double));
    gblcblk   = (double*)realloc(gblcblk,  n_layers*8*sizeof(double));
  }
  MMM2D_self_energy();
}


void  MMM2D_dielectric_layers_force_contribution()
{
  int c, i, j;
  Cell  *celll;
  int      npl;
  Particle *pl, *p1;
  double dist2, d[3];
  double charge_factor;
  double a[3];
  double force[3]={0, 0, 0};

  if (!mmm2d_params.dielectric_contrast_on) return;

  if(this_node==0) {
    c=1;
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;
        
    for(i = 0; i < npl; i++) {
      // printf("enter mmm2d_dielectric %f \n",my_left[2]);
      force[0]=0;force[1]=0;force[2]=0;
      p1 = &pl[i];
      for(j = 0; j < npl; j++) {
	a[0]=pl[j].r.p[0]; a[1]=pl[j].r.p[1];a[2]=-pl[j].r.p[2];
       	layered_get_mi_vector(d, p1->r.p, a);
       	dist2 = sqrlen(d);
	charge_factor=p1->p.q*pl[j].p.q*mmm2d_params.delta_mid_bot;
	add_mmm2d_coulomb_pair_force(charge_factor, d, sqrt(dist2), dist2, force);
      }
      for (j = 0; j < 3; j++) { 
	p1->f.f[j] += force[j];
      }
    }
  }
  
  if(this_node==n_nodes-1) {
  
    c=n_layers;
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;
    for(i = 0; i < npl; i++) {
      force[0]=0;force[1]=0;force[2]=0;
      p1 = &pl[i];
      for(j = 0; j < npl; j++) {
	a[0]=pl[j].r.p[0];  a[1]=pl[j].r.p[1]; a[2]=2*box_l[2]-pl[j].r.p[2];
	layered_get_mi_vector(d, p1->r.p, a);
	dist2 = sqrlen(d);
	charge_factor=p1->p.q*pl[j].p.q*mmm2d_params.delta_mid_top; 
	add_mmm2d_coulomb_pair_force(charge_factor, d, sqrt(dist2), dist2, force);
      }
      for (j = 0; j < 3; j++) { 
	p1->f.f[j] += force[j];
      }
    }
  }

}

double  MMM2D_dielectric_layers_energy_contribution()
{
  int c, i, j;
  Cell  *celll;
  int      npl;
  Particle *pl, *p1;
  double dist2, d[3];
  double charge_factor;
  double a[3];
  double eng=0.0;
  // prefactor for the charged plate interaction removal
  double corr_pref = coulomb.prefactor*C_2PI*ux*uy;

  if (!mmm2d_params.dielectric_contrast_on) return 0.0;

  if(this_node==0) {
    c=1;
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;
    
    for(i = 0; i < npl; i++) {
      p1 = &pl[i];
      for(j = 0; j < npl; j++) {
	a[0]=pl[j].r.p[0];  a[1]=pl[j].r.p[1];a[2]=-pl[j].r.p[2];
       	layered_get_mi_vector(d, p1->r.p, a);
       	dist2 = sqrlen(d);
	charge_factor = mmm2d_params.delta_mid_bot*p1->p.q*pl[j].p.q;
	eng+=mmm2d_coulomb_pair_energy(charge_factor, d, dist2, sqrt(dist2)) + corr_pref*charge_factor*d[2];
      }
    }
  }
  
  if(this_node==n_nodes-1) {
    c=n_layers;
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;
    for(i = 0; i < npl; i++) {
      p1 = &pl[i];
      for(j = 0; j < npl; j++) {
	a[0]=pl[j].r.p[0];  a[1]=pl[j].r.p[1]; a[2]=2*box_l[2]-pl[j].r.p[2];
	layered_get_mi_vector(d, p1->r.p, a);
	dist2 = sqrlen(d);
	charge_factor=mmm2d_params.delta_mid_top*p1->p.q*pl[j].p.q;
	eng+=mmm2d_coulomb_pair_energy(charge_factor, d, dist2, sqrt(dist2)) - corr_pref*charge_factor*d[2];
      }
    }
  }
  return 0.5*eng;
}

#endif
