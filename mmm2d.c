#include <math.h>
#include <mpi.h>
#include "communication.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "cells.h"
#include "cells.h"
#include "config.h"
#include "mmm2d.h"
#include "mmm-common.h"
#include "utils.h"
#include "specfunc.h"
#include "layered.h"

#ifdef ELECTROSTATICS

// #define CHECKPOINTS
#if 0
#define LOG_FORCES(x) x
#else
#define LOG_FORCES(x)
#endif

/****************************************
 * LOCAL DEFINES
 ****************************************/

/* Largest reasonable cutoff for far formula */
#define MAXIMAL_FAR_CUT 50

/* Largest reasonable cutoff for Bessel function */
#define MAXIMAL_B_CUT 20

/* Largest reasonable order of polygamma series */
#define MAXIMAL_POLYGAMMA 50

/* number of steps in the complex cutoff table */
#define COMPLEX_STEP 16
/* map numbers from 0 to 1/sqrt(2) onto the complex cutoff table
   (with security margin) */
#define COMPLEX_FAC (COMPLEX_STEP/(1/M_SQRT2 + 0.01))

/****************************************
 * LOCAL VARIABLES
 ****************************************/

/* up to that error the sums in the NF are evaluated */
static double part_error;

/* cutoffs for the bessel sum */
static IntList besselCutoff = {NULL, 0, 0};

/* cutoffs for the complex sum */
static int  complexCutoff[COMPLEX_STEP];
static DoubleList  bon = {NULL, 0, 0};

static double ux, ux2, uy, uy2, uz;
static double cell_h;
static double self_energy;

MMM2D_struct mmm2d_params = { 1e100, 10, 1 };

/****************************************
 * LOCAL ARRAYS
 ****************************************/

/* product decomposition data organization. For the cell blocks
   it is assumed that the lower blocks part is in the lower half.
   This has to have positive sign, so that has to be first. */
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

/* number of local particles */
static int n_localpart = 0;

/* temporary buffers for product decomposition */
static double *partblk = NULL;
/* for all local cells including ghosts */
static double *lclcblk = NULL;
/* collected data from the cells above the top neighbor
   of a cell rsp. below the bottom neighbor
   (P=below, M=above, as the signs in the exp). */
static double *gblcblk = NULL;

/* structure for storing of sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/* sin/cos caching */ 
static SCCache *scxcache = NULL;
static int    n_scxcache;  
static SCCache *scycache = NULL;
static int    n_scycache;  

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

/* NEAR FORMULA */
/****************/

/* complex evaluation */
static void prepareBernoulliNumbers(int nmax);

/* cutoff error setup. Returns error message or NULL */
static char *MMM2D_tune_near(double error);

/** energy of all local particles with their copies */
void MMM2D_self_energy();

/* FAR FORMULA */
/***************/
/* sin/cos storage */
static void prepare_scx_cache();
static void prepare_scy_cache();
/* common code */
static void distribute(int size);
/* 2 pi sign(z) code */
static void add_2pi_signz();
static double twopi_z_energy();
/* p=0 per frequency code */
static void setup_P(int p, double omega);
static void add_P_force();
static double   P_energy(double omega);
/* q=0 per frequency code */
static void setup_Q(int q, double omega);
static void add_Q_force();
static double   Q_energy(double omega);
/* p,q <> 0 per frequency code */
static void setup_PQ(int p, int q, double omega);
static void add_PQ_force(int p, int q, double omega);
static double   PQ_energy(double omega);

/* cutoff error setup. Returns error message or NULL */
static char *MMM2D_tune_far(double error);

/* COMMON */
/**********/

void MMM2D_setup_constants()
{
  ux  = 1/box_l[0];
  ux2 = ux*ux;
  uy  = 1/box_l[1];
  uy2 = uy*uy;  
  uz  = 1/box_l[2];

  switch (cell_structure.type) {
  case CELL_STRUCTURE_NSQUARE:
    cell_h = box_l[2];
    break;
  case CELL_STRUCTURE_LAYERED:
    cell_h = layer_h;
    break;
  default:
    cell_h = 0;
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

MDINLINE void clear_vec(double *pdc, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc[i] = 0;
}

MDINLINE void copy_vec(double *pdc_d, double *pdc_s, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s[i];
}

MDINLINE void add_vec(double *pdc_d, double *pdc_s1, double *pdc_s2, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc_d[i] = pdc_s1[i] + pdc_s2[i];
}

MDINLINE void scale_vec(double scale, double *pdc, int size)
{
  int i;
  for (i = 0; i < size; i++)
    pdc[i] *= scale;
}

MDINLINE double *block(double *p, int index, int size)
{
  return &p[index*size];
}

MDINLINE double *blwentry(double *p, int index, int e_size)
{
  return &p[2*index*e_size];
}

MDINLINE double *abventry(double *p, int index, int e_size)
{
  return &p[(2*index + 1)*e_size];
}

void distribute(int e_size)
{
  int c, node, inv_node;
  double sendbuf[8];
  double recvbuf[8];
  MPI_Status status;

  if (this_node == 0) {
    /* no neighbor below */
    clear_vec(blwentry(gblcblk, 0, e_size), e_size);
    clear_vec(blwentry(lclcblk, 0, e_size), e_size);
  }
  if (this_node == n_nodes - 1) {
    /* no neighbor above */
    clear_vec(abventry(gblcblk, n_layers - 1, e_size), e_size);
    clear_vec(abventry(lclcblk, n_layers + 1, e_size), e_size);
  }

  /* send/recv to/from other nodes. Also builds up the gblcblk. */
  for (node = 0; node < n_nodes; node++) {
    inv_node = n_nodes - node - 1;
    /* up */
    if (node == this_node) {
      /* calculate sums of cells below */
      for (c = 1; c < n_layers; c++)
	add_vec(blwentry(gblcblk, c, e_size), blwentry(gblcblk, c - 1, e_size), blwentry(lclcblk, c - 1, e_size), e_size);
 
      /* calculate my ghost contribution only if a node below exists */
      if (node + 1 < n_nodes) {
	add_vec(sendbuf, blwentry(gblcblk, n_layers - 1, e_size), blwentry(lclcblk, n_layers - 1, e_size), e_size);
	copy_vec(sendbuf + e_size, blwentry(lclcblk, n_layers, e_size), e_size);
	MPI_Send(sendbuf, 2*e_size, MPI_DOUBLE, node + 1, 0, MPI_COMM_WORLD);
      }
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*e_size, MPI_DOUBLE, node, 0, MPI_COMM_WORLD, &status);
      copy_vec(blwentry(gblcblk, 0, e_size), recvbuf, e_size);
      copy_vec(blwentry(lclcblk, 0, e_size), recvbuf + e_size, e_size);
    }

    /* up */
    if (inv_node == this_node) {
      /* calculate sums of all cells above */
      for (c = n_layers + 1; c > 2; c--)
	add_vec(abventry(gblcblk, c - 3, e_size), abventry(gblcblk, c - 2, e_size), abventry(lclcblk, c, e_size), e_size);

      /* calculate my ghost contribution only if a node above exists */
      if (inv_node -  1 >= 0) {
	add_vec(sendbuf, abventry(gblcblk, 0, e_size), abventry(lclcblk, 2, e_size), e_size);
	copy_vec(sendbuf + e_size, abventry(lclcblk, 1, e_size), e_size);
	MPI_Send(sendbuf, 2*e_size, MPI_DOUBLE, inv_node - 1, 0, MPI_COMM_WORLD);
      }
    }
    else if (inv_node - 1 == this_node) {
      MPI_Recv(recvbuf, 2*e_size, MPI_DOUBLE, inv_node, 0, MPI_COMM_WORLD, &status);
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

static void add_2pi_signz()
{
  int np, c, i;
  double pref = coulomb.prefactor*C_2PI*ux*uy, add;
  Particle *part;
  int e_size = 1, size = 2;
  double *othcblk;

  /* calculate local cellblks. partblks don't make sense */
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    lclcblk[size*c] = 0;
    for (i = 0; i < np; i++)
      lclcblk[size*c] += part[i].p.q;
    lclcblk[size*c] *= pref;
    /* just to be able to use the standard distribution. Here below
       and above terms are the same */
    lclcblk[size*c + 1] = lclcblk[size*c];
  }

  distribute(e_size);

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

static double twopi_z_energy()
{
  int np, c, i;
  double pref = -coulomb.prefactor*C_2PI*ux*uy;
  Particle *part;
  int e_size = 2, size = 4;
  double *othcblk;
  double eng = 0;

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

  distribute(e_size);

  for (c = 1; c <= n_layers; c++) {
    othcblk = block(gblcblk, c - 1, size);
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++)
      eng += part[i].p.q*(part[i].r.p[2]*othcblk[ABEQQP] - othcblk[ABEQZP] -
			  part[i].r.p[2]*othcblk[ABEQQM] + othcblk[ABEQZM]);
  }
  return eng;
}

/*****************************************************************/
/* PoQ exp sum */
/*****************************************************************/

static void setup_P(int p, double omega)
{
  int np, c, i, ic, o = (p-1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*4*M_PI*ux*uy;
  double e;
  double *llclcblk;
  int size = 4;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*part[i].r.p[2]);

      partblk[size*ic + POQESM] = part[i].p.q*scxcache[o + ic].s/e;
      partblk[size*ic + POQESP] = part[i].p.q*scxcache[o + ic].s*e;
      partblk[size*ic + POQECM] = part[i].p.q*scxcache[o + ic].c/e;
      partblk[size*ic + POQECP] = part[i].p.q*scxcache[o + ic].c*e;

      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, llclcblk, size);
  }
}

static void setup_Q(int q, double omega)
{
  int np, c, i, ic, o = (q-1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*4*M_PI*ux*uy;
  double e;
  double *llclcblk;
  int size = 4;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*part[i].r.p[2]);

      partblk[size*ic + POQESM] = part[i].p.q*scycache[o + ic].s/e;
      partblk[size*ic + POQESP] = part[i].p.q*scycache[o + ic].s*e;
      partblk[size*ic + POQECM] = part[i].p.q*scycache[o + ic].c/e;
      partblk[size*ic + POQECP] = part[i].p.q*scycache[o + ic].c*e;

      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, llclcblk, size);
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
  Particle *part;
  double *othcblk;
  int size = 4;
  double eng = 0;
  double pref = 1/omega;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
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
  Particle *part;
  double *othcblk;
  int size = 4;
  double eng = 0;
  double pref = 1/omega;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
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

static void setup_PQ(int p, int q, double omega)
{
  int np, c, i, ic, ox = (p - 1)*n_localpart, oy = (q - 1)*n_localpart;
  Particle *part;
  double pref = coulomb.prefactor*8*M_PI*ux*uy;
  double e;
  double *llclcblk;
  int size = 8;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    llclcblk = block(lclcblk, c, size);

    clear_vec(llclcblk, size);

    for (i = 0; i < np; i++) {
      e = exp(omega*part[i].r.p[2]);
      
      partblk[size*ic + PQESSM] = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q/e;
      partblk[size*ic + PQESCM] = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q/e;
      partblk[size*ic + PQECSM] = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q/e;
      partblk[size*ic + PQECCM] = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q/e;

      partblk[size*ic + PQESSP] = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e;
      partblk[size*ic + PQESCP] = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e;
      partblk[size*ic + PQECSP] = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e;
      partblk[size*ic + PQECCP] = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e;

      add_vec(llclcblk, llclcblk, block(partblk, ic, size), size);
      ic++;
    }
    scale_vec(pref, llclcblk, size);
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
  int np, c, i, ic;
  Particle *part;
  double *othcblk;
  int size = 8;
  double eng = 0;
  double pref = 1/omega;

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    othcblk = block(gblcblk, c - 1, size);

    for (i = 0; i < np; i++) {
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

void MMM2D_add_far_force()
{
  int p, q;
  double omega;

  if (cell_structure.type != CELL_STRUCTURE_LAYERED || n_nodes*n_layers < 3)
    return;

  add_2pi_signz();
  checkpoint("************2piz", 0, 0, 1);

  prepare_scx_cache();
  prepare_scy_cache();

  /* the second condition is just for the case of numerical accident */
  for (p = 1; ux*(p - 1) < mmm2d_params.far_cut && p <= n_scxcache; p++) {
    omega = C_2PI*ux*p;
    setup_P(p, omega);
    distribute(2);
    add_P_force();
    checkpoint("************distri p", p, 0, 2);
  }

  for (q = 1; uy*(q - 1) < mmm2d_params.far_cut && q <= n_scycache; q++) {
    omega = C_2PI*uy*q;
    setup_Q(q, omega);
    distribute(2);
    add_Q_force();
    checkpoint("************distri q", 0, q, 2);
  }

  for (p = 1; p - 1 < ux*mmm2d_params.far_cut  && p <= n_scxcache ; p++) {
    for (q = 1; SQR(ux*(p - 1)) + SQR(uy*(q - 1)) < mmm2d_params.far_cut2 && q <= n_scycache; q++) {
      omega = C_2PI*sqrt(SQR(ux*p) + SQR(uy*q));
      setup_PQ(p, q, omega);
      distribute(4);
      add_PQ_force(p, q, omega);
      checkpoint("************distri pq", p, q, 4);
    }
  }
}

double MMM2D_far_energy()
{
  double eng;
  int p, q;
  double omega;

  // It's not really far...
  eng = self_energy;

  if (cell_structure.type != CELL_STRUCTURE_LAYERED || n_nodes*n_layers < 3)
    return 0.5*eng;

  eng += twopi_z_energy();
  checkpoint("E************2piz", 0, 0, 2);

  prepare_scx_cache();
  prepare_scy_cache();

  /* the second condition is just for the case of numerical accident */
  for (p = 1; ux*(p - 1) < mmm2d_params.far_cut && p <= n_scxcache; p++) {
    omega = C_2PI*ux*p;
    setup_P(p, omega);
    distribute(2);
    eng += P_energy(omega);
    checkpoint("E************distri p", p, 0, 2);
  }

  for (q = 1; uy*(q - 1) < mmm2d_params.far_cut && q <= n_scycache; q++) {
    omega = C_2PI*uy*q;
    setup_Q(q, omega);
    distribute(2);
    eng += Q_energy(omega);
    checkpoint("E************distri q", 0, q, 2);
  }

  for (p = 1; p - 1 < ux*mmm2d_params.far_cut  && p <= n_scxcache ; p++) {
    for (q = 1; SQR(ux*(p - 1)) + SQR(uy*(q - 1)) < mmm2d_params.far_cut2 && q <= n_scycache; q++) {
      omega = C_2PI*sqrt(SQR(ux*p) + SQR(uy*q));
      setup_PQ(p, q, omega);
      distribute(4);
      eng += PQ_energy(omega);
      checkpoint("E************distri pq", p, q, 4);
    }
  }

  /* we count both i<->j and j<->i, so return just half of it */
  return 0.5*eng;
}

static char *MMM2D_tune_far(double error)
{
  double err, step;
  mmm2d_params.far_cut = 2;
  /* so it is at least at one border optimal */
  step = dmin(ux, uy);
  do {
    err = exp(-2*M_PI*mmm2d_params.far_cut*cell_h)/cell_h*
      (C_2PI*mmm2d_params.far_cut + 2*(ux + uy) + 1/cell_h);
    mmm2d_params.far_cut += step;
  }
  while (err > error && mmm2d_params.far_cut < MAXIMAL_FAR_CUT);
  if (mmm2d_params.far_cut >= MAXIMAL_FAR_CUT)
    return "Far cutoff too large, decrease n_layers or the error bound";
  // fprintf(stderr, "far cutoff %g %g\n", mmm2d_params.far_cut, err);
  mmm2d_params.far_cut -= step;
  mmm2d_params.far_cut2 = SQR(mmm2d_params.far_cut);
  return NULL;
}

/****************************************
 * NEAR FORMULA
 ****************************************/

static char *MMM2D_tune_near(double error)
{
  int P, n, i;
  double uxrho2m2max, uxrhomax2;
  int p;
  double T, pref, err, exponent;
  double L, sum;

  /* yes, it's y only... */
  if (cell_h > box_l[1]/2)
    return "Layer height too large for MMM2D near formula, increase n_layers";
  if (ux*box_l[1] >= 3/M_SQRT2 )
    return "box_l[1]/box_l[0] too large for MMM2D near formula, please exchange x and y";

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
    return "Could find not reasonable Bessel cutoff. Please decrease n_layers or the error bound";
  // fprintf(stderr, "bessel cutoff %d %g\n", P, err);

  realloc_intlist(&besselCutoff, besselCutoff.n = P);
  for (p = 1; p < P; p++)
    besselCutoff.e[p-1] = (int)floor(((double)P)/(2*p)) + 1;

  /* complex sum, determine cutoffs (dist dependent) */
  T = log(part_error/(16*M_SQRT2)*box_l[1]*box_l[2]);
  for (i = 0; i < COMPLEX_STEP; i++)
    complexCutoff[i] = (int)ceil(T/log((i+1)/COMPLEX_FAC));
  prepareBernoulliNumbers(complexCutoff[COMPLEX_STEP-1]);

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
    return "Could find not reasonable Polygamma cutoff. Consider exchanging x and y";
  // fprintf(stderr, "polygamma cutoff %d %g\n", n, err);

  return NULL;
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
      bon.e[l-1] =  2.0/l;
    else
      bon.e[l-1] = -2.0/l;      
  }
}

void add_mmm2d_coulomb_pair_force(Particle *p1, Particle *p2,
				  double d[3], double dl2, double dl)
{
  double F[3];
  double pref = coulomb.prefactor*p1->p.q*p2->p.q;
  double z2   = d[2]*d[2];
  double rho2 = d[1]*d[1] + z2;
  int i;

  if (pref != 0.0) {
    F[0] = F[1] = F[2] = 0;

    /* Bessel sum */
    {
      int p, l;
      double k1;
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
	  k0Sum  += K0(freq*rho_l);
	  k1 = K1(freq*rho_l)/rho_l;
	  k1Sum  += k1;
	  k1ySum += k1*ypl;

	  ypl   = d[1] - l*box_l[1];
	  rho_l = sqrt(ypl*ypl + z2);
	  k0Sum  += K0(freq*rho_l);
	  k1 = K1(freq*rho_l)/rho_l;
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

      end = complexCutoff[(int)ceil(COMPLEX_FAC*uy2*rho2)];
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

    for (i = 0; i < 3; i++) {
      p1->f.f[i] += pref*F[i];
      p2->f.f[i] -= pref*F[i];
    }
  }
}

MDINLINE double calc_mmm2d_copy_pair_energy(double d[3])
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

    end = complexCutoff[(int)ceil(COMPLEX_FAC*uy2*rho2)];
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

double mmm2d_coulomb_pair_energy(Particle *p1, Particle *p2,
				 double dv[3], double d2, double d)
{
  double eng, pref = coulomb.prefactor*p1->p.q*p2->p.q;

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

  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;

    for (i = 0; i < np; i++) {
      self_energy += seng*SQR(part[i].p.q);
    }
  }
}

/****************************************
 * COMMON PARTS
 ****************************************/

void MMM2D_init()
{
  char *err;

  if (!PERIODIC(0) || !PERIODIC(1) || PERIODIC(2)) {
    fprintf(stderr, "ERROR: MMM2D requires periodicity 1 1 0\n");
    errexit();
  }
  
  if (cell_structure.type != CELL_STRUCTURE_LAYERED &&
      cell_structure.type != CELL_STRUCTURE_NSQUARE) {
    fprintf(stderr, "ERROR: MMM2D requires layered (or n-square) cellsystem\n");
    errexit();
  }

  MMM2D_setup_constants();
  if ((err = MMM2D_tune_near(mmm2d_params.maxPWerror))) {
    fprintf(stderr, "ERROR: MMM2D auto-retuning: %s\n", err);
    errexit();
  }
  if (mmm2d_params.far_calculated) {
    if ((err = MMM2D_tune_far(mmm2d_params.maxPWerror))) {
      fprintf(stderr, "ERROR: MMM2D auto-retuning: %s\n", err);
      errexit();
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
    scxcache = realloc(scxcache, n_scxcache*n_localpart*sizeof(SCCache));
    scycache = realloc(scycache, n_scycache*n_localpart*sizeof(SCCache));
    
    partblk   = realloc(partblk,  n_localpart*8*sizeof(double));
    lclcblk   = realloc(lclcblk,  n_cells*8*sizeof(double));
    gblcblk   = realloc(gblcblk,  n_layers*8*sizeof(double));
  }
  MMM2D_self_energy();
}

int set_mmm2d_params(Tcl_Interp *interp, double maxPWerror, double far_cut)
{
  char *err;

  mmm2d_params.maxPWerror = maxPWerror;

  MMM2D_setup_constants();

  if ((err = MMM2D_tune_near(maxPWerror))) {
    coulomb.method = COULOMB_NONE;
    Tcl_AppendResult(interp, err, NULL);
    return TCL_ERROR;
  }

  if (far_cut >= 0) {
    mmm2d_params.far_cut = far_cut;
    mmm2d_params.far_cut2 = SQR(far_cut);
    mmm2d_params.far_calculated = 0;
  }
  else {
    if ((err = MMM2D_tune_far(maxPWerror))) {
      coulomb.method = COULOMB_NONE;
      Tcl_AppendResult(interp, err, NULL);
      return TCL_ERROR;
    }
    mmm2d_params.far_calculated = 1;    
  }

  coulomb.method = COULOMB_MMM2D;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

#endif
