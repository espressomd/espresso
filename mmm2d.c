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

/****************************************
 * LOCAL DEFINES
 ****************************************/

/* Largest reasonable cutoff for Bessel function */
#define MAXIMAL_B_CUT 30

/* number of steps in the complex cutoff table */
#define COMPLEX_STEP 16
/* map numbers from 0 to 0.5 onto the complex cutoff table
   (with security margin) */
#define COMPLEX_FAC (COMPLEX_STEP/0.501)

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
static double cell_h, cell_h_i;

MMM2D_struct mmm2d_params = { 1e100, 10, 1 };

/****************************************
 * LOCAL ARRAYS
 ****************************************/

/* structure for storing of sin and cos values */
typedef struct {
  double s, c;
} SCCache;

/* structure for product decomposition data.
   letters:
   PQ:  1: sin/cos x
        2: sin/cos y
*/
typedef struct {
  double ss, sc;
  double cs, cc;
} PQCacheE;

/* full cache with plus amd minus sign of exp for PQ */
typedef struct {
  PQCacheE p;
  PQCacheE m;
} PQCache;

/* number of local particles */
static int n_localpart = 0;

/* temporary buffers for product decomposition */
static PQCache *pqpartblk = NULL;
/* for all local cells including ghosts */
static PQCache *pqlclcblk = NULL;
/* collected data from the cells above (m) rsp. below (p). m and p refers to the
   exp sign, which has to be positive for cells below. */
static PQCache  pqothcblk;
/* collected data from the cells above. To avoid allocating sabvcblk too
   often, declared here, but used only in add_*_blocks. */
PQCacheE *pqsavcblk = NULL;

/* cache for the P or Q loops */
typedef SCCache PoQCacheE;

/* full cache with plus amd minus sign of exp for P or Q */
typedef struct {
  PoQCacheE p;
  PoQCacheE m;
} PoQCache;

/* temporary buffers for product decomposition */
static PoQCache *poqpartblk = NULL;
/* for all local cells including ghosts */
static PoQCache *poqlclcblk = NULL;
/* collected data from the cells above (m) rsp. below (p). m and p refers to the
   exp sign, which has to be positive for cells below. */
static PoQCache  poqothcblk;
/* collected data from the cells above. To avoid allocating sabvcblk too
   often, declared here, but used only in add_*_blocks. */
PoQCacheE *poqsavcblk = NULL;

/* sin/cos caching */ 
static SCCache *scxcache = NULL;
static SCCache *scycache = NULL;

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

/* NEAR FORMULA */
/****************/
/* Bessel evaluation */
static void prepareBesselCutoffs(int P);

/* complex evaluation */
static void prepareBernoulliNumbers(int nmax);

/* cutoff error setup */
static void MMM2D_tune_near(double error);

/* FAR FORMULA */
/***************/
/* sin/cos storage */
static void prepare_scx_cache();
static void prepare_scy_cache();
/* 2 pi sign(z) code */
static void add_2pi_signz();
/* common code */
static void sum_and_distribute_PoQ(double scale);
static void sum_and_distribute_PQ(double scale);
/* p=0 per frequency code */
static void setup_P(int p, double omega);
static void add_P_force(double scale);
/* q=0 per frequency code */
static void setup_Q(int q, double omega);
static void add_Q_force(double scale);
/* p,q <> 0 per frequency code */
static void setup_PQ(int p, int q, double omega);
static void add_PQ_force(double scale, int p, int q, double omega);

/* cutoff error setup */
static void MMM2D_tune_far(double error);

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
  cell_h_i = 1/cell_h;
}

/****************************************
 * FAR FORMULA
 ****************************************/

static void prepare_scx_cache()
{
  int np, c, i, ic, freq, o, max = (int)ceil(ux*mmm2d_params.far_cut) + 1;
  double pref, arg;
  Particle *part;
  
  for (freq = 1; freq <= max; freq++) {
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
  int np, c, i, ic, freq, o, max = (int)ceil(uy*mmm2d_params.far_cut) + 1;
  double pref, arg;
  Particle *part;
  
  for (freq = 1; freq <= max; freq++) {
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
/* 2 pi sign(z) */
/*****************************************************************/

static void add_2pi_signz()
{
  int np, c, i;
  int node, inv_node;
  int cnt_below, cnt_above;
  int sendbuf;
  double add;
  Particle *part;
  MPI_Status status;

  if (this_node == 0)
    /* now neighbor below */
    cnt_below = 0;
  if (this_node == n_nodes - 1)
    /* now neighbor above */
    cnt_above = 0;

  /* send/recv to/from other nodes, no periodicity */
  for (node = 0; node < n_nodes - 1; node++) {
    inv_node = n_nodes - node - 1;
    /* up */
    if (node == this_node) {
      /* calculate sums of all cells below */
      sendbuf = cnt_below;
      /* all my cells except top one */
      for (c = 1; c < n_layers; c++)
	sendbuf += cells[c].n;
      /* no need to send the top cell count, this is already
	 done by the ghost mechanism */
      MPI_Send(&sendbuf, 1, MPI_INT, node + 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node)
      MPI_Recv(&cnt_below, 1, MPI_INT, node - 1, 0, MPI_COMM_WORLD, &status);

    /* down */
    if (node == this_node) {
      /* calculate sums of all cells above */
      sendbuf = cnt_above;
      /* all my cells except bottom one */
      for (c = 2; c <= n_layers; c++)
	sendbuf += cells[c].n;
      MPI_Send(&sendbuf, 1, MPI_INT, node - 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node)
      MPI_Recv(&cnt_above, 1, MPI_INT, node - 1, 0, MPI_COMM_WORLD, &status);
  }

  /* run through from bottom layer. This means we first add up all top contributions and
     subtract them later again */
  for (c = n_layers + 1; c > 2; c--)
    cnt_above += cells[c].n;

  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    add = C_2PI*ux*uy*(cnt_above - cnt_below);
    for (i = 0; i < np; i++)
      part[i].f.f[2] += add;

    // fprintf(stderr, "cell %d(%d part) ca %d cb %d\n", c, np, cnt_above, cnt_below);

    if (c < n_layers) {
      cnt_below += cells[c - 1].n;
      cnt_above -= cells[c + 2].n;
    }
  }
}

/*****************************************************************/
/* PoQ exp sum */
/*****************************************************************/

static void setup_P(int p, double omega)
{
  int np, c, i, ic, o = (p-1)*n_localpart;
  Particle *part;
  double e;
  double box_b;

  box_b = 0;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - box_b));

      poqpartblk[ic].m.s = part[i].p.q*scxcache[o + ic].s/e;
      poqpartblk[ic].p.s = part[i].p.q*scxcache[o + ic].s*e;
      poqpartblk[ic].m.c = part[i].p.q*scxcache[o + ic].c/e;
      poqpartblk[ic].p.c = part[i].p.q*scxcache[o + ic].c*e;

      ic++;
    }
    box_b += layer_h;
  }
}

static void setup_Q(int q, double omega)
{
  int np, c, i, ic, o = (q-1)*n_localpart;
  Particle *part;
  double e;
  double box_b;

  box_b = 0;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - box_b));

      poqpartblk[ic].m.s = part[i].p.q*scycache[o + ic].s/e;
      poqpartblk[ic].p.s = part[i].p.q*scycache[o + ic].s*e;
      poqpartblk[ic].m.c = part[i].p.q*scycache[o + ic].s/e;
      poqpartblk[ic].p.c = part[i].p.q*scycache[o + ic].s*e;

      ic++;
    }
    box_b += layer_h;
  }
}

MDINLINE void clear_PoQCacheE(PoQCacheE *pdc)
{
  pdc->s = 0;
  pdc->c = 0;
}

MDINLINE void copy_PoQCacheE(PoQCacheE *pdc_d, PoQCacheE *pdc_s)
{
  pdc_d->s = pdc_s->s;
  pdc_d->c = pdc_s->c;
}

MDINLINE void copy_scaled_PoQCacheE(PoQCacheE *pdc_d, double scale, PoQCacheE *pdc_s)
{
  pdc_d->s = pdc_s->s*scale;
  pdc_d->c = pdc_s->c*scale;
}

MDINLINE void add_PoQCacheE(PoQCacheE *pdc_d, PoQCacheE *pdc_s)
{
  pdc_d->s += pdc_s->s;
  pdc_d->c += pdc_s->c;
}

MDINLINE void scale_and_add_PoQCacheE(double scale, PoQCacheE *pdc_d, PoQCacheE *pdc_s)
{
  pdc_d->s = pdc_d->s*scale + pdc_s->s;
  pdc_d->c = pdc_d->c*scale + pdc_s->c;
}

MDINLINE void scale_and_add_scaled_PoQCacheE(double scale_d, PoQCacheE *pdc_d, double scale_s, PoQCacheE *pdc_s)
{
  pdc_d->s = pdc_d->s*scale_d + pdc_s->s*scale_s;
  pdc_d->c = pdc_d->c*scale_d + pdc_s->c*scale_s;
}

MDINLINE void copyadd_scaled_PoQCacheE(PoQCacheE *pdc_d, double scale_1, PoQCacheE *pdc_s1, double scale_2, PoQCacheE *pdc_s2)
{
  pdc_d->s = pdc_s1->s*scale_1 + pdc_s2->s*scale_2;
  pdc_d->c = pdc_s1->c*scale_1 + pdc_s2->c*scale_2;
}

static void sum_and_distribute_PoQ(double scale)
{
  int np, c, i, ic;
  int node, inv_node;
  PoQCacheE recvbuf[2];
  PoQCacheE sendbuf[2];
  MPI_Status status;

  /* calculate local cellblks */
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    clear_PoQCacheE(&poqlclcblk[c].p);
    clear_PoQCacheE(&poqlclcblk[c].m);
    for (i = 0; i < np; i++) {
      add_PoQCacheE(&poqlclcblk[c].m, &poqpartblk[ic].m);
      add_PoQCacheE(&poqlclcblk[c].p, &poqpartblk[ic].p);

      ic++;
    }
  }

  if (this_node == 0) {
    /* no neighbor below */
    clear_PoQCacheE(&poqothcblk.p);
    clear_PoQCacheE(&poqlclcblk[0].p);
  }
  if (this_node == n_nodes - 1) {
    /* no neighbor above */
    clear_PoQCacheE(&poqothcblk.m);
    clear_PoQCacheE(&poqlclcblk[n_layers + 1].m);
  }

  /* send/recv to/from other nodes, no periodicity */
  for (node = 0; node < n_nodes - 1; node++) {
    inv_node = n_nodes - node - 1;
    /* up, plus sign */
    if (node == this_node) {
      /* calculate sums of all cells below, p sign */
      copy_PoQCacheE(&sendbuf[0], &poqothcblk.p);
      /* all my cells except top one */
      for (c = 1; c < n_layers; c++)
	scale_and_add_PoQCacheE(scale, &sendbuf[0], &poqlclcblk[c].p);

      copy_PoQCacheE(&sendbuf[1], &poqlclcblk[n_layers].p);

      MPI_Send(sendbuf, 2*sizeof(PoQCacheE), MPI_BYTE, node + 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*sizeof(PoQCacheE), MPI_BYTE, node - 1, 0, MPI_COMM_WORLD, &status);
      copy_PoQCacheE(&poqothcblk.p, &recvbuf[0]);
      copy_PoQCacheE(&poqlclcblk[0].p, &recvbuf[1]);
    }

    /* down, minus sign */
    if (node == this_node) {
      /* calculate sums of all cells above, m sign */
      copy_PoQCacheE(&sendbuf[0], &poqothcblk.m);
      /* all my cells except bottom one */
      for (c = 2; c <= n_layers; c++)
	scale_and_add_PoQCacheE(scale, &sendbuf[0], &poqlclcblk[c].m);

      copy_PoQCacheE(&sendbuf[1], &poqlclcblk[1].m);

      MPI_Send(sendbuf, 2*sizeof(PoQCacheE), MPI_BYTE, node - 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*sizeof(PoQCacheE), MPI_BYTE, node + 1, 0, MPI_COMM_WORLD, &status);
      copy_PoQCacheE(&poqothcblk.m, &recvbuf[0]);
      copy_PoQCacheE(&poqlclcblk[n_layers + 1].m, &recvbuf[1]);
    }
  }
}

static void add_P_force(double scale)
{
  int np, c, i, ic;
  Particle *part;
  double pref_scale2 = 4*M_PI*ux*uy*scale*scale;
  PoQCacheE blwcblk, *avcblk;

  /* calculate the above blocks, relative to my layers bottoms. */
  copy_scaled_PoQCacheE(&poqsavcblk[n_layers - 1], pref_scale2, &poqothcblk.m);
  for (c = n_layers - 1; c >= 1; c--)
    copyadd_scaled_PoQCacheE(&poqsavcblk[c - 1], scale, &poqsavcblk[c], pref_scale2, &poqlclcblk[c + 2].m);
  /* blwcblk is done on the fly */
  copy_scaled_PoQCacheE(&blwcblk, pref_scale2, &poqothcblk.p);

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    avcblk = &poqsavcblk[c - 1];

    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      part[i].f.f[0] +=
	poqpartblk[ic].m.s*blwcblk.c - poqpartblk[ic].m.c*blwcblk.s +
	poqpartblk[ic].p.s*avcblk->c - poqpartblk[ic].p.c*avcblk->s;
      part[i].f.f[2] +=
	poqpartblk[ic].m.c*blwcblk.c + poqpartblk[ic].m.s*blwcblk.s -
	poqpartblk[ic].p.c*avcblk->c - poqpartblk[ic].p.s*avcblk->s;

      ic++;
    }

    /* update blwcblk */
    scale_and_add_scaled_PoQCacheE(scale, &blwcblk, pref_scale2, &poqlclcblk[c - 1].p);
  }
}

static void add_Q_force(double scale)
{
  int np, c, i, ic;
  Particle *part;
  double pref_scale2 = 4*M_PI*ux*uy*scale*scale;
  PoQCacheE blwcblk, *avcblk;

  /* calculate the above blocks, relative to my layers bottoms. */
  copy_scaled_PoQCacheE(&poqsavcblk[n_layers - 1], pref_scale2, &poqothcblk.m);
  for (c = n_layers - 1; c >= 1; c--)
    copyadd_scaled_PoQCacheE(&poqsavcblk[c - 1], scale, &poqsavcblk[c], pref_scale2, &poqlclcblk[c + 2].m);
  /* blwcblk is done on the fly */
  copy_scaled_PoQCacheE(&blwcblk, pref_scale2, &poqothcblk.p);

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    avcblk = &poqsavcblk[c - 1];

    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      part[i].f.f[1] +=
	poqpartblk[ic].m.s*blwcblk.c - poqpartblk[ic].m.c*blwcblk.s +
	poqpartblk[ic].p.s*avcblk->c - poqpartblk[ic].p.c*avcblk->s;
      part[i].f.f[2] +=
	poqpartblk[ic].m.c*blwcblk.c + poqpartblk[ic].m.s*blwcblk.s -
	poqpartblk[ic].p.c*avcblk->c - poqpartblk[ic].p.s*avcblk->s;

      ic++;
    }

    /* update blwcblk */
    scale_and_add_scaled_PoQCacheE(scale, &blwcblk, pref_scale2, &poqlclcblk[c - 1].p);
  }
}

/*****************************************************************/
/* PQ particle blocks */
/*****************************************************************/

static void setup_PQ(int p, int q, double omega)
{
  int np, c, i, ic, ox = (p - 1)*n_localpart, oy = (q - 1)*n_localpart;
  Particle *part;
  double e;
  double box_b;

  box_b = 0;
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      e = exp(omega*(part[i].r.p[2] - box_b));
      
      pqpartblk[ic].m.ss = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q/e;
      pqpartblk[ic].m.sc = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q/e;
      pqpartblk[ic].m.cs = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q/e;
      pqpartblk[ic].m.cc = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q/e;

      pqpartblk[ic].p.ss = scxcache[ox + ic].s*scycache[oy + ic].s*part[i].p.q*e;
      pqpartblk[ic].p.sc = scxcache[ox + ic].s*scycache[oy + ic].c*part[i].p.q*e;
      pqpartblk[ic].p.cs = scxcache[ox + ic].c*scycache[oy + ic].s*part[i].p.q*e;
      pqpartblk[ic].p.cc = scxcache[ox + ic].c*scycache[oy + ic].c*part[i].p.q*e;

      ic++;
    }
    box_b += layer_h;
  }
}

MDINLINE void clear_PQCacheE(PQCacheE *pdc)
{
  pdc->ss = pdc->sc =
    pdc->cs = pdc->cc = 0;
}

MDINLINE void copy_PQCacheE(PQCacheE *pdc_d, PQCacheE *pdc_s)
{
  pdc_d->ss = pdc_s->ss;
  pdc_d->sc = pdc_s->sc;
  pdc_d->cs = pdc_s->cs;
  pdc_d->cc = pdc_s->cc;
}

MDINLINE void copy_scaled_PQCacheE(PQCacheE *pdc_d, double scale, PQCacheE *pdc_s)
{
  pdc_d->ss = pdc_s->ss*scale;
  pdc_d->sc = pdc_s->sc*scale;
  pdc_d->cs = pdc_s->cs*scale;
  pdc_d->cc = pdc_s->cc*scale;
}

MDINLINE void add_PQCacheE(PQCacheE *pdc_d, PQCacheE *pdc_s)
{
  pdc_d->ss += pdc_s->ss;
  pdc_d->sc += pdc_s->sc;
  pdc_d->cs += pdc_s->cs;
  pdc_d->cc += pdc_s->cc;
}

MDINLINE void scale_and_add_PQCacheE(double scale, PQCacheE *pdc_d, PQCacheE *pdc_s)
{
  pdc_d->ss = pdc_d->ss*scale + pdc_s->ss;
  pdc_d->sc = pdc_d->sc*scale + pdc_s->sc;
  pdc_d->cs = pdc_d->cs*scale + pdc_s->cs;
  pdc_d->cc = pdc_d->cc*scale + pdc_s->cc;
}

MDINLINE void scale_and_add_scaled_PQCacheE(double scale_d, PQCacheE *pdc_d, double scale_s, PQCacheE *pdc_s)
{
  pdc_d->ss = pdc_d->ss*scale_d + pdc_s->ss*scale_s;
  pdc_d->sc = pdc_d->sc*scale_d + pdc_s->sc*scale_s;
  pdc_d->cs = pdc_d->cs*scale_d + pdc_s->cs*scale_s;
  pdc_d->cc = pdc_d->cc*scale_d + pdc_s->cc*scale_s;
}

MDINLINE void copyadd_scaled_PQCacheE(PQCacheE *pdc_d, double scale_1, PQCacheE *pdc_s1, double scale_2, PQCacheE *pdc_s2)
{
  pdc_d->ss = pdc_s1->ss*scale_1 + pdc_s2->ss*scale_2;
  pdc_d->sc = pdc_s1->sc*scale_1 + pdc_s2->sc*scale_2;
  pdc_d->cs = pdc_s1->cs*scale_1 + pdc_s2->cs*scale_2;
  pdc_d->cc = pdc_s1->cc*scale_1 + pdc_s2->cc*scale_2;
}

static void sum_and_distribute_PQ(double scale)
{
  int np, c, i, ic;
  int node, inv_node;
  PQCacheE recvbuf[2];
  PQCacheE sendbuf[2];
  MPI_Status status;

  /* calculate local cellblks */
  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    np   = cells[c].n;
    clear_PQCacheE(&pqlclcblk[c].p);
    clear_PQCacheE(&pqlclcblk[c].m);
    for (i = 0; i < np; i++) {
      add_PQCacheE(&pqlclcblk[c].m, &pqpartblk[ic].m);
      add_PQCacheE(&pqlclcblk[c].p, &pqpartblk[ic].p);

      ic++;
    }
  }

  if (this_node == 0) {
    /* no neighbor below */
    clear_PQCacheE(&pqothcblk.p);
    clear_PQCacheE(&pqlclcblk[0].p);
  }
  if (this_node == n_nodes - 1) { 
    /* no neighbor above */
    clear_PQCacheE(&pqothcblk.m);
    clear_PQCacheE(&pqlclcblk[n_layers + 1].m);
  }

  /* send/recv to/from other nodes, no periodicity */
  for (node = 0; node < n_nodes - 1; node++) {
    inv_node = n_nodes - node - 1;
    /* up, plus sign */
    if (node == this_node) {
      /* calculate sums of all cells below, p sign */
      copy_PQCacheE(&sendbuf[0], &pqothcblk.p);
      /* all my cells except top one */
      for (c = 1; c < n_layers; c++)
	scale_and_add_PQCacheE(scale, &sendbuf[0], &pqlclcblk[c].p);

      copy_PQCacheE(&sendbuf[1], &pqlclcblk[n_layers].p);

      MPI_Send(sendbuf, 2*sizeof(PQCacheE), MPI_BYTE, node + 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*sizeof(PQCacheE), MPI_BYTE, node - 1, 0, MPI_COMM_WORLD, &status);
      copy_PQCacheE(&pqothcblk.p, &recvbuf[0]);
      copy_PQCacheE(&pqlclcblk[0].p, &recvbuf[1]);
    }

    /* down, minus sign */
    if (node == this_node) {
      /* calculate sums of all cells above, m sign */
      copy_PQCacheE(&sendbuf[0], &pqothcblk.m);
      /* all my cells except bottom one */
      for (c = 2; c <= n_layers; c++)
	scale_and_add_PQCacheE(scale, &sendbuf[0], &pqlclcblk[c].m);

      copy_PQCacheE(&sendbuf[1], &pqlclcblk[1].m);

      MPI_Send(sendbuf, 2*sizeof(PQCacheE), MPI_BYTE, node - 1, 0, MPI_COMM_WORLD);
    }
    else if (node + 1 == this_node) {
      MPI_Recv(recvbuf, 2*sizeof(PQCacheE), MPI_BYTE, node + 1, 0, MPI_COMM_WORLD, &status);
      copy_PQCacheE(&pqothcblk.m, &recvbuf[0]);
      copy_PQCacheE(&pqlclcblk[n_layers + 1].m, &recvbuf[1]);
    }
  }
}

static void add_PQ_force(double scale, int p, int q, double omega)
{
  int np, c, i, ic;
  Particle *part;
  double pref_scale2 = 16*M_PI*M_PI*ux*uy*scale*scale;
  double pref_x = ux*p/omega;
  double pref_y = uy*q/omega;
  PQCacheE blwcblk, *avcblk;

  /* calculate the above blocks, relative to my layers bottoms. */
  copy_scaled_PQCacheE(&pqsavcblk[n_layers - 1], pref_scale2, &pqothcblk.m);
  for (c = n_layers - 1; c >= 1; c--)
    copyadd_scaled_PQCacheE(&pqsavcblk[c - 1], scale, &pqsavcblk[c], pref_scale2, &pqlclcblk[c + 2].m);
  /* blwcblk is done on the fly */
  copy_scaled_PQCacheE(&blwcblk, pref_scale2, &pqothcblk.p);

  ic = 0;
  for (c = 1; c <= n_layers; c++) {
    avcblk = &pqsavcblk[c - 1];

    np   = cells[c].n;
    part = cells[c].part;
    for (i = 0; i < np; i++) {
      part[i].f.f[0] +=
	pref_x*(pqpartblk[ic].m.sc*blwcblk.cc + pqpartblk[ic].m.ss*blwcblk.cs -
		pqpartblk[ic].m.cc*blwcblk.sc - pqpartblk[ic].m.cs*blwcblk.ss +
		pqpartblk[ic].p.sc*avcblk->cc + pqpartblk[ic].p.ss*avcblk->cs -
		pqpartblk[ic].p.cc*avcblk->sc - pqpartblk[ic].p.cs*avcblk->ss);
      part[i].f.f[1] +=
	pref_y*(pqpartblk[ic].m.cs*blwcblk.cc - pqpartblk[ic].m.cc*blwcblk.cs +
		pqpartblk[ic].m.ss*blwcblk.sc - pqpartblk[ic].m.sc*blwcblk.ss +
		pqpartblk[ic].p.cs*avcblk->cc - pqpartblk[ic].p.cc*avcblk->cs +
		pqpartblk[ic].p.ss*avcblk->sc - pqpartblk[ic].p.sc*avcblk->ss);
      part[i].f.f[2] +=
	       (pqpartblk[ic].m.cc*blwcblk.cc + pqpartblk[ic].m.cs*blwcblk.cs +
	        pqpartblk[ic].m.sc*blwcblk.sc + pqpartblk[ic].m.ss*blwcblk.ss -
	        pqpartblk[ic].p.cc*avcblk->cc - pqpartblk[ic].p.cs*avcblk->cs -
	        pqpartblk[ic].p.sc*avcblk->sc - pqpartblk[ic].p.ss*avcblk->ss);

      ic++;
    }

    /* update blwcblk */
    scale_and_add_scaled_PQCacheE(scale, &blwcblk, pref_scale2, &pqlclcblk[c - 1].p);
  }
}

void MMM2D_add_far()
{
  int p, q;
  double scale, omega;

  if (cell_structure.type != CELL_STRUCTURE_LAYERED || n_layers < 3)
    return;

  add_2pi_signz();

  prepare_scx_cache();
  prepare_scy_cache();

  for (p = 1; p - 1 < ux*mmm2d_params.far_cut; p++) {
    omega = C_2PI*ux*p; 
    scale = exp(-omega*cell_h);
    setup_P(p, omega);
    sum_and_distribute_PoQ(scale);
    add_P_force(scale);
  }

  for (q = 1; q - 1 < uy*mmm2d_params.far_cut; q++) {
    omega = C_2PI*uy*q; 
    scale = exp(-omega*cell_h);
    setup_Q(p, omega);
    sum_and_distribute_PoQ(scale);
    add_Q_force(scale);
  }

  for (p = 1; p - 1 < ux*mmm2d_params.far_cut; p++) {
    for (q = 1; SQR(ux*(p - 1)) + SQR(uy*(q - 1)) < SQR(mmm2d_params.far_cut); q++) {
      omega = C_2PI*sqrt(SQR(ux*p) + SQR(uy*q)); 
      scale = exp(-omega*cell_h);
      setup_PQ(p, q, omega);
      sum_and_distribute_PQ(scale);
      add_PQ_force(scale, p, q, omega);
    }
  }
 
}

static void MMM2D_tune_far(double error)
{
  //  fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!! fixed far cutoff ten\n");
  mmm2d_params.far_cut = 10;
  mmm2d_params.far_calculated = 1;
}

/****************************************
 * NEAR FORMULA
 ****************************************/

static void MMM2D_tune_near(double error)
{
  int P, n, i;
  double T, err;
  double uxrho2m2max, uxrhomax2;

  mmm2d_params.maxPWerror = error;

  /* yes, it's y only... */
  if (cell_h > box_l[1]/2) {
    fprintf(stderr, "layer height too large for MMM2D near formula, increase n_layers\n");
    errexit();
  }
  if (ux*box_l[1]/M_SQRT2 > 1.5) {
    fprintf(stderr, "box_l[1]/box_l[0] too large for MMM2D near formula\n");
    errexit();
  }

  /* error is split into three parts:
     one part for bessel, one for complex
     and one for polygamma cutoff */
  part_error = error/3;

  /* Bessel sum, determine cutoff */
  P = 1;
  T  = exp(M_PI)/M_PI;
  do {
    int p;
    double L, sum;
    L = P*M_PI;
    sum = 0;
    for (p = 1; p < P; p++)
      sum += p*exp(-M_PI*p);
    err = 16*M_PI*exp(-L)*(T*((L + 1)/M_PI - 1) + sum);
    P++;
  }
  while (err > part_error);
  P--;
  prepareBesselCutoffs(P);
  
  /* complex sum, determine cutoffs (dist dependent) */
  T = log(part_error/16/M_SQRT2*box_l[1]*box_l[2]);
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
  while (err > part_error);
}

static void prepareBesselCutoffs(int P)
{
  int p;
  double L2MPI;

  realloc_intlist(&besselCutoff, besselCutoff.n = P);
  L2MPI = 0.5*P;
  for (p = 1; p < P; p++)
    besselCutoff.e[p-1] = (int)ceil(L2MPI/p) + 1;
}

static void prepareBernoulliNumbers(int bon_order)
{
  int l;
  double pref, fac;

  /* Bernoulli over faculty */
  static double bon_table[15] = {
    1.00000000000000000000000000000,
    0.083333333333333333333333333333,
    -0.00138888888888888888888888888889,
    0.000033068783068783068783068783069,
    -8.2671957671957671957671957672e-07,
    2.0876756987868098979210090321e-08,
    -5.2841901386874931848476822022e-10,
    1.3382536530684678832826980975e-11,
    -3.3896802963225828668301953912e-13,
    8.5860620562778445641359054504e-15,
    -2.1748686985580618730415164239e-16,
    5.5090028283602295152026526089e-18,
    -1.3954464685812523340707686264e-19,
    3.5347070396294674716932299778e-21,
    -8.9535174270375468504026113181e-23
  };

  if (bon_order < 2)
    bon_order = 2;

  realloc_doublelist(&bon, bon.n = bon_order);
  
  fac  = C_2PI*C_2PI*uy2;
  /* the ux is multiplied in to bessel, complex and psi at once, not here */
  pref = 2*fac;
  for(l = 1; (l <= bon_order) && (l < 15); l++) {
    bon.e[l-1] = pref*bon_table[l];
    pref *= fac;
  }
  for (; l <= bon_order; l++) {
    if (l & 1)
      bon.e[l-1] =  2.0/l;
    else
      bon.e[l-1] = -2.0/l;      
  }
}

MDINLINE void calc_mmm2d_copy_pair_force(double d[3], double *F)
{
  double z2     = d[2]*d[2];
  double rho2   = d[1]*d[1] + z2;
  int i;

  /* Bessel sum */
  {
    int p, l;
    double k1;
    double k0Sum, k1ySum, k1Sum;
    double freq;
    double rho_l, ypl;
    double c, s;

    F[0] = F[1] = F[2] = 0;

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

    ztn_r = zeta_r = d[2];
    ztn_i = zeta_i = d[1];
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
  }
  // fprintf(stderr, "explcit force %f %f %f\n", F[0], F[1], F[2]);
}

void add_mmm2d_coulomb_pair_force(Particle *p1, Particle *p2,
				  double dv[3], double d2, double d)
{
  int i;
  double ri3;
  double F[3], pref = p1->p.q*p2->p.q;

  if (pref != 0.0) {
    calc_mmm2d_copy_pair_force(dv, F);

    ri3 = 1/(d2*d);
    for (i = 0; i < 3; i++) {
      p1->f.f[i] += pref*(F[i] + ri3*dv[i]);
      p2->f.f[i] -= pref*(F[i] + ri3*dv[i]);
    }
    // fprintf(stderr, "%d<->%d: %f %f %f\n", p1->p.identity, p2->p.identity, dv[0], dv[1], dv[2]);
  }
}

/****************************************
 * COMMON PARTS
 ****************************************/

void MMM2DdoInternalChecks()
{
  FILE *out;
  int j;
  double Wn, jfac;

  out = fopen("bon.log", "w");
  fprintf(stderr, "generating bon.log for %d bernoulli numbers\n",
	  bon.n);
  Wn = C_2PI*C_2PI;
  jfac = 2;
  for (j = 1; j <= bon.n; j++) {
    fprintf(out, "%d %f %f\n", j, bon.e[j-1], jfac*bon.e[j-1]/(2*Wn)*(2*j));
    Wn *= C_2PI*C_2PI;
    jfac *= (2*j+1)*(2*j+2);
  }
  fclose(out);

  out = fopen("bessel.log", "w");
  fprintf(stderr, "generating bessel.log for a cutoff of %d\n", besselCutoff.n);
  for (j = 1; j < besselCutoff.n; j++)
    fprintf(out, "%d %d\n", j, besselCutoff.e[j - 1]);
  fclose(out);

  /*
  fprintf(stderr, "psi order is %d\n", psiMaxOrder);
  for (n = 0; n < psiMaxOrder; n++) {
    char buf[256];
    sprintf(buf, "psi%02d.log", n);
    fprintf(stderr, "generating %s for order %d\n", buf, n);
    out = fopen(buf, "w");
    printf("generating %s\n", buf);
    for (i = 0; i <= psiEven[n].order; i++)
      fprintf(out, "%d %f\n", i, psiEven[n].coeff[i]);

    for (i = 0; i <= psiOdd[n].order; i++)
      fprintf(out, "%d %f\n", i, psiOdd[n].coeff[i]);
    fclose(out);
  }
  */
}

void MMM2D_init()
{
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
  MMM2D_tune_near(mmm2d_params.maxPWerror);
  if (mmm2d_params.far_calculated)
    MMM2D_tune_far(mmm2d_params.maxPWerror);

  /* allocate caches */
  n_localpart = cells_get_n_particles();
  scxcache = realloc(scxcache, (int)(ceil(ux*mmm2d_params.far_cut) + 1)*n_localpart*sizeof(SCCache));
  scycache = realloc(scycache, (int)(ceil(uy*mmm2d_params.far_cut) + 1)*n_localpart*sizeof(SCCache));

  pqpartblk   = realloc(pqpartblk,  n_localpart*sizeof(PQCache));
  pqlclcblk   = realloc(pqlclcblk,  n_cells*sizeof(PQCache));
  pqsavcblk   = realloc(pqsavcblk,  n_layers*sizeof(PQCacheE));
  
  poqpartblk  = realloc(poqpartblk, n_localpart*sizeof(PoQCache));
  poqlclcblk  = realloc(poqlclcblk, n_cells*sizeof(PoQCache));
  poqsavcblk  = realloc(poqsavcblk, n_layers*sizeof(PoQCacheE));
}

int set_mmm2d_params(Tcl_Interp *interp, double maxPWerror)
{
  MMM2D_setup_constants();
  MMM2D_tune_near(maxPWerror);
  MMM2D_tune_far(maxPWerror);
  coulomb.method = COULOMB_MMM2D;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

#endif
