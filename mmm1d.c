#include "mmm1d.h"
#include "polynom.h"
#include "specfunc.h"
#include "communication.h"
#include "cells.h"
#include "grid.h"
#include "tuning.h"
#include "utils.h"
#include <mpi.h>
#include <tcl.h>

#ifdef ELECTROSTATICS


/*
  cwz-build-command: gmake
odd strategy:

    0 1 2 3 4
  0 0 1 0 3 0
  1 1 1 2 1 4
  2 0 2 2 3 2
  3 3 1 3 3 4
  4 0 4 2 4 4

    0 1 2 3 4
0.
s   0 1 2 3 4
r   0 1 2 3 4
1. offset=1
s   1 2 3 4 0
r   4 0 1 2 3
2. offset=3
s   3 4 0 1 2
r   2 3 4 0 1

  0 1 2 3 4
0 0 1 2 2 1
1   0 1 2 2
2     0 1 2
3       0 1
4         0

even strategy:

    0 1 2 3
  0 0 1 0 3
  1 1 1 2 1
  2 0 2 2 3
  3 3 1 3 3

    0 1 2 3
0.
s   0 1 2 3
r   0 1 2 3
1. offset = 1
s   1 2 3 *
r   * 0 1 2
2. offset = 3
s   3 * 0 1
r   2 3 * 0

  0 1 2 3
0 0 1 2 2
1   0 1 2
2     0 1
3       0

*/

#define C_2PI     (2*M_PI)
#define C_GAMMA   0.57721566490153286060651209008
#define C_2LOG4PI -5.0620484939385815859557831885
#define C_2PISQR  C_2PI*C_2PI
/* precision of polygamma functions. More is unnecessary, the Bessel
   functions are not better anyways... */
#define EPS 1e-10

/* How many trial calculations */
#define TEST_INTEGRATIONS 1000

/* Largest reasonable cutoff for Bessel function */
#define MAXIMAL_B_CUT 30

/* Granularity of the radius scan in multiples of box_l[2] */
#define RAD_STEPPING 0.1

/** modified polygamma functions. See Arnold,Holm 2002 */
static Polynom *modPsi = NULL;
static int      n_modPsi = 0, eo_n_nodes = 0;
static double L_i, L2, L2_i, prefL2_i, prefL3_i;

MMM1D_struct mmm1d_params = { 0.05, 5, 1e-5, 1, 1};

MDINLINE double mod_psi_even(int n, double x)
{ return evaluateAsTaylorSeriesAt(&modPsi[2*n],x*x); }

MDINLINE double mod_psi_odd(int n, double x)
{ return x*evaluateAsTaylorSeriesAt(&modPsi[2*n+1], x*x); }

static void preparePolygammaEven(int n, double binom, Polynom *series)
{
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv = 2*n;
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    maxx = 0.25;
    alloc_doublelist(series, 1);
    series->e[0] = 2*(1 - C_GAMMA);
    for (order = 1;; order += 1) {
      x_order = 2*order;
      coeff = -2*hzeta(x_order + 1, 2);
      if (fabs(maxx*coeff)*(4.0/3.0) < EPS)
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = coeff;
      maxx *= 0.25;
    }
    series->n = order;
  }
  else {
    // even, n > 0
    maxx = 1;
    pref = 2;
    init_doublelist(series);
    for (order = 0;; order++) {
      // only even exponents of x
      x_order = 2*order;
      coeff = pref*hzeta(1 + deriv + x_order, 2);
      if ((fabs(maxx*coeff)*(4.0/3.0) < EPS) && (x_order > deriv))
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = -binom*coeff;
      maxx *= 0.25;
      pref *= (1.0 + deriv/(x_order + 1));
      pref *= (1.0 + deriv/(x_order + 2));
    }
    series->n = order;
  }
}

static void preparePolygammaOdd(int n, double binom, Polynom *series)
{
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv  = 2*n + 1;
  maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  pref = 2*deriv*(1 + deriv);
  init_doublelist(series);
  for (order = 0;; order++) {
    // only odd exponents of x
    x_order = 2*order + 1;
    coeff = pref*hzeta(1 + deriv + x_order, 2);
    if ((fabs(maxx*coeff)*(4.0/3.0) < EPS) && (x_order > deriv))
      break;
    realloc_doublelist(series, order + 1);
    series->e[order] = -binom*coeff;
    maxx *= 0.25;
    pref *= (1.0 + deriv/(x_order + 1));
    pref *= (1.0 + deriv/(x_order + 2));
  }
  series->n = order;
}

double determine_bessel_cutoff(double switch_rad, double maxPWerror, int maxP)
{
  /* this calculates an upper bound to all force components and the potential */

  // fprintf(stderr, "determ: %f %f %d\n", switch_rad, maxPWerror, maxP);

  double err;
  double Lz = box_l[2];
  double fac = 2*M_PI/Lz*switch_rad;
  double pref = 4/Lz*dmax(1, 2*M_PI/Lz);
  int P = (int)ceil(3*box_l[2]/2/M_PI/switch_rad);
  do {
    err = pref*exp(-fac*(P-1))*(1 + P*(exp(fac) - 1))/SQR(1 - exp(fac));
    // fprintf(stderr, "%d %e\n", P, err); */
    P++;
  } while (err > maxPWerror && P <= maxP);
  P--;
  return P;
}

int set_mmm1d_params(Tcl_Interp *interp, double bjerrum, double switch_rad,
		     int bessel_cutoff, double maxPWerror)
{
  char buffer[32 + 2*TCL_DOUBLE_SPACE];
  double int_time, min_time=1e200, min_rad = -1;
  double maxrad = box_l[2];
  /* more doesn't work with the current implementation as N_psi = 2 fixed */

  if (bessel_cutoff < 0 && switch_rad < 0) {
    /* determine besselcutoff and optimal switching radius */
    mmm1d_params.bjerrum = 1;
    for (switch_rad = RAD_STEPPING*box_l[2]; switch_rad < maxrad; switch_rad += RAD_STEPPING*box_l[2]) {
      mmm1d_params.bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
      /* no reasonable cutoff possible */
      if (mmm1d_params.bessel_cutoff == MAXIMAL_B_CUT)
	continue;
      mmm1d_params.far_switch_radius_2 = switch_rad*switch_rad;

      /* initialize mmm1d temporary structures */
      mpi_bcast_coulomb_params();
      /* perform force calculation test */
      int_time = time_force_calc(TEST_INTEGRATIONS);

      sprintf(buffer, "r= %f c= %d t= %f ms\n", switch_rad, mmm1d_params.bessel_cutoff, int_time);
      // fprintf(stderr, buffer);
      Tcl_AppendResult(interp, buffer, (char*)NULL);

      if (int_time < min_time) {
	min_time = int_time;
	min_rad = switch_rad;
      }
      /* stop if all hope is vain... */
      else if (int_time > 2*min_time)
	break;
    }
    if (min_rad < 0) {
      fprintf(stderr, "set_mmm1d_params: internal error");
      errexit();
    }
    switch_rad    = min_rad;
    bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
  }

  if (bessel_cutoff < 0) {
    /* determine besselcutoff to achieve at least the given pairwise error */
    bessel_cutoff = determine_bessel_cutoff(switch_rad, maxPWerror, MAXIMAL_B_CUT);
    if (bessel_cutoff == MAXIMAL_B_CUT) {
      Tcl_AppendResult(interp, "could not find reasonable bessel cutoff", (char *)NULL);
      return TCL_ERROR;
    }
  }

  if (switch_rad <= 0 || switch_rad > box_l[2]) {
    Tcl_AppendResult(interp, "switching radius is not between 0 and box_l[2]", (char *)NULL);
    return TCL_ERROR;
  }
  if (bessel_cutoff <=0) {
    Tcl_AppendResult(interp, "bessel cutoff too small", (char *)NULL);
    return TCL_ERROR;
  }

  mmm1d_params.bjerrum = bjerrum;
  mmm1d_params.far_switch_radius_2 = switch_rad*switch_rad;
  mmm1d_params.bessel_cutoff = bessel_cutoff;

  mmm1d_params.maxPWerror = maxPWerror;

  mpi_bcast_coulomb_params();

  return TCL_OK;
}

void MMM1D_recalcTables()
{
  /* polygamma, determine order */
  int n;
  double binom, err;
  double rho2m2max;

  if (modPsi != NULL) {
    for (n = 0; n < 2*n_modPsi; n++)
      realloc_doublelist(&modPsi[n], 0);
    modPsi = NULL;
    n_modPsi = 0;
  }

  n = 0;
  binom = 1.0;
  rho2m2max = 1.0;
  do {
    n_modPsi++;
    modPsi = realloc(modPsi, 2*n_modPsi*sizeof(Polynom));
    
    preparePolygammaEven(n, binom, &modPsi[2*n_modPsi - 2]);
    preparePolygammaOdd(n, binom, &modPsi[2*n_modPsi - 1]);
    
    err = fabs(2*mod_psi_even(n,rho2m2max));
    rho2m2max *= 0.5;
    binom *= (-0.5 - n)/(double)(n+1);
    n++;
  }
  while (err > 0.1*mmm1d_params.maxPWerror);
}

void MMM1D_init()
{
  /* round up n_nodes to next odd number */
  eo_n_nodes = (n_nodes/2)*2 + 1;

  /* precalculate some constants */
  L_i  = 1/box_l[2];
  L2   = box_l[2]*box_l[2];
  L2_i = L_i*L_i;
  prefL2_i = mmm1d_params.prefactor*L2_i;
  prefL3_i = prefL2_i*L_i;

  MMM1D_recalcTables();
}

MDINLINE void calc_pw_force(double dx, double dy, double dz,
			    double *Fx, double *Fy, double *Fz)
{
  double rxy2, rxy2_d, z_d;
  double r2, r, pref;

  dz -= rint(dz/box_l[2])*box_l[2];
  
  rxy2   = dx*dx + dy*dy;
  rxy2_d = rxy2*L2_i;
  z_d    = dz*L_i;
  r2     = rxy2 + dz*dz;
  r      = sqrt(r2);

  if (rxy2 <= mmm1d_params.far_switch_radius_2) {
    /* near range formula */
    double sr, sz, r2nm1, rt, rt2, shift_z;
    int n;

    /* polygamma summation */
    sr = 0;
    sz = mod_psi_odd(0, z_d);

    r2nm1 = 1.0;
    for (n = 1; n < n_modPsi; n++) {
      double deriv = 2*n;
      double mpe   = mod_psi_even(n, z_d);
      double mpo   = mod_psi_odd(n, z_d);
      double r2n   = r2nm1*rxy2_d;

      sz +=         r2n*mpo;
      sr += deriv*r2nm1*mpe;

      if (fabs(deriv*r2nm1*mpe) < mmm1d_params.maxPWerror)
	break;

      r2nm1 = r2n;
    }

    *Fx = prefL3_i*sr*dx;
    *Fy = prefL3_i*sr*dy;
    *Fz = prefL2_i*sz;

    /* real space parts */

    pref = mmm1d_params.prefactor/(r2*r); 
    *Fx += pref*dx;
    *Fy += pref*dy;
    *Fz += pref*dz;

    shift_z = dz + box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    rt  = sqrt(rt2);
    pref = mmm1d_params.prefactor/(rt2*rt); 
    *Fx += pref*dx;
    *Fy += pref*dy;
    *Fz += pref*shift_z;

    shift_z = dz - box_l[2];
    rt2 = rxy2 + shift_z*shift_z;
    rt  = sqrt(rt2);
    pref = mmm1d_params.prefactor/(rt2*rt); 
    *Fx += pref*dx;
    *Fy += pref*dy;
    *Fz += pref*shift_z;
  }
  else {
    /* far range formula */
    double rxy   = sqrt(rxy2);
    double rxy_d = rxy*L_i;
    double sr = 0, sz = 0;
    int bp;

    for (bp = 1; bp < mmm1d_params.bessel_cutoff; bp++) {
      double fq = C_2PI*bp;    
      sr += bp*K1(fq*rxy_d)*cos(fq*z_d);
      sz += bp*K0(fq*rxy_d)*sin(fq*z_d);
    }
    sr *= L2_i*4*C_2PI;
    sz *= L2_i*4*C_2PI;
    
    pref = mmm1d_params.prefactor*(sr/rxy + 2*L_i/rxy2);

    *Fx = pref*dx;
    *Fy = pref*dy;
    *Fz = mmm1d_params.prefactor*sz;
  }
}

void MMM1D_calc_forces()
{
  MPI_Status status;

  double *send_coords;
  double *send_forces;
  int n_send_coords;
  int send_node;

  double *recv_coords = NULL;
  double *recv_forces = NULL;
  int n_recv_coords = 0, max_recv_coords = 0;
  int recv_node;

  int m1, n1, o1, p1;
  int m2, n2, o2, p2;
  Particle *part1, *part2, *p_part1, *p_part2;
  int n_part1, n_part2;
  int ind1, ind2;

  int *sizes;
  int buffer_part;
  int offset;
  double chpref;
  double dx, dy, dz;
  double Fx, Fy, Fz;

  if (mmm1d_params.maxPWerror <= 0)
    return;

  if (periodic[0] != 0 || periodic[1] != 0 || periodic[2] != 1) {
    fprintf(stderr, "MMM1D must have periodicity 0,0,1 !\n");
    errexit();
  }

  /* prepare send and fetch buffers */
  n_send_coords = cells_get_n_particles();
  send_coords = malloc(n_send_coords*sizeof(double)*4);
  send_forces = malloc(n_send_coords*sizeof(double)*3);
  sizes = malloc(sizeof(int)*n_nodes);
  MPI_Allgather(&n_send_coords, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

  /* interaction on node itself and setup of send buffer */
  buffer_part = 0;
  INNER_CELLS_LOOP(m1, n1, o1) {
    ind1 = CELL_IND(m1,n1,o1);
    n_part1 = cells[ind1].pList.n;
    p_part1 = cells[ind1].pList.part;
    for (p1 = 0; p1 < n_part1; p1++) {
      part1 = &p_part1[p1];

      /* append particle to send buffer */
      send_coords[4*buffer_part    ] = part1->r.p[0];
      send_coords[4*buffer_part + 1] = part1->r.p[1];
      send_coords[4*buffer_part + 2] = part1->r.p[2];
      send_coords[4*buffer_part + 3] = part1->r.q;
      /*
	fprintf(stderr, "%d: send %f %f %f %f\n", this_node,
	send_coords[4*buffer_part    ],
	send_coords[4*buffer_part + 1],
	send_coords[4*buffer_part + 2],
	send_coords[4*buffer_part + 3]);
      */
      buffer_part++;

      /* interactions with all the other particles on this node */
      INNER_CELLS_LOOP(m2, n2, o2) {
	ind2 = CELL_IND(m2,n2,o2);
	if (ind2 < ind1)
	  continue;
	n_part2 = cells[ind2].pList.n;
	p_part2 = cells[ind2].pList.part;
	for (p2 = 0; p2 < n_part2; p2++) {
	  part2 = &p_part2[p2];

	  /* no self energy */
	  if (ind1 == ind2 && part2->r.identity <= part1->r.identity)
	    continue;

	  chpref = part1->r.q*part2->r.q;

	  if (chpref != 0.0) {
	    dx = part1->r.p[0] - part2->r.p[0];
	    dy = part1->r.p[1] - part2->r.p[1];
	    dz = part1->r.p[2] - part2->r.p[2];

	    calc_pw_force(dx, dy, dz, &Fx, &Fy, &Fz);

	    Fx *= chpref;
	    Fy *= chpref;
	    Fz *= chpref;

	    part1->f[0] += Fx;
	    part1->f[1] += Fy;
	    part1->f[2] += Fz;
	    part2->f[0] -= Fx;
	    part2->f[1] -= Fy;
	    part2->f[2] -= Fz;
	  }
	}
      }
    }
  }

  for (offset = 1; offset < eo_n_nodes; offset += 2) {
    send_node = (this_node + offset) % eo_n_nodes;
    recv_node = (this_node - offset + eo_n_nodes) % eo_n_nodes;
    if (recv_node < n_nodes)
      n_recv_coords = sizes[recv_node];
    else
      n_recv_coords = 0;
    /*
      fprintf(stderr, "round %2d: %2d(%2d) sn %2d rn %2d(%2d)\n",
      offset, this_node,
      n_send_coords, send_node,
      recv_node, n_recv_coords);
    */

    if (n_recv_coords > max_recv_coords) {
      max_recv_coords = n_recv_coords;
      recv_coords = realloc(recv_coords, max_recv_coords*sizeof(double)*4);
      recv_forces = realloc(recv_forces, max_recv_coords*sizeof(double)*3);
    }
    /* send away my particles coordinates and charge */
    if (send_node < n_nodes)
      MPI_Send(send_coords, 4*n_send_coords, MPI_DOUBLE, send_node,
	       1, MPI_COMM_WORLD);
    /* and get the particles which are my duty now, and calculate */
    if (recv_node < n_nodes) {
      MPI_Recv(recv_coords, 4*n_recv_coords, MPI_DOUBLE, recv_node,
	       1, MPI_COMM_WORLD, &status);
      for (p2 = 0; p2 < 3*n_recv_coords; p2++)
	recv_forces[p2] = 0;

      INNER_CELLS_LOOP(m1, n1, o1) {
	ind1 = CELL_IND(m1,n1,o1);
	n_part1 = cells[ind1].pList.n;
	p_part1 = cells[ind1].pList.part;
	for (p1 = 0; p1 < n_part1; p1++) {
	  part1 = &p_part1[p1];
	  for (p2 = 0; p2 < n_recv_coords; p2++) {
	    /*
	      fprintf(stderr, "%d: got %f %f %f %f\n", this_node,
	      recv_coords[4*p2    ],
	      recv_coords[4*p2 + 1],
	      recv_coords[4*p2 + 2],
	      recv_coords[4*p2 + 3]);
	    */
	    dx = part1->r.p[0] - recv_coords[4*p2    ];
	    dy = part1->r.p[1] - recv_coords[4*p2 + 1];
	    dz = part1->r.p[2] - recv_coords[4*p2 + 2];
	    chpref = part1->r.q*recv_coords[4*p2 + 3];

	    if (chpref != 0.0) {
	      calc_pw_force(dx, dy, dz, &Fx, &Fy, &Fz);

	      Fx *= chpref;
	      Fy *= chpref;
	      Fz *= chpref;
	      part1->f[0] += Fx;
	      part1->f[1] += Fy;
	      part1->f[2] += Fz;
	      recv_forces[3*p2    ] -= Fx;
	      recv_forces[3*p2 + 1] -= Fy;
	      recv_forces[3*p2 + 2] -= Fz;
	    }
	    /*
	      fprintf(stderr, "%d: sendf %f %f %f\n", this_node,
	      recv_forces[3*p2    ],
	      recv_forces[3*p2 + 1],
	      recv_forces[3*p2 + 2]);
	    */
	  }
	}
      }

      /* send away other nodes forces */
      if (recv_node < n_nodes)
	MPI_Send(recv_forces, 3*n_recv_coords, MPI_DOUBLE, recv_node,
		 1, MPI_COMM_WORLD);
    }

    /* and collect my ones */
    if (send_node < n_nodes) {
      MPI_Recv(send_forces, 3*n_send_coords, MPI_DOUBLE, send_node,
	       1, MPI_COMM_WORLD, &status);    
      buffer_part = 0;
      INNER_CELLS_LOOP(m1, n1, o1) {
	ind1 = CELL_IND(m1,n1,o1);
	n_part1 = cells[ind1].pList.n;
	p_part1 = cells[ind1].pList.part;
	for (p1 = 0; p1 < n_part1; p1++) {
	  part1 = &p_part1[p1];
	  part1->f[0] += send_forces[3*buffer_part    ];
	  part1->f[1] += send_forces[3*buffer_part + 1];
	  part1->f[2] += send_forces[3*buffer_part + 2];
	  /*
	    fprintf(stderr, "%d: gotf %f %f %f\n", this_node,
	    send_forces[3*buffer_part    ],
	    send_forces[3*buffer_part + 1],
	    send_forces[3*buffer_part + 2]);
	  */
	  buffer_part++;
	}
      }
    }
  }

  /*
  INNER_CELLS_LOOP(m1, n1, o1) {
    ind1 = CELL_IND(m1,n1,o1);
    n_part1 = cells[ind1].pList.n;
    p_part1 = cells[ind1].pList.part;
    for (p1 = 0; p1 < n_part1; p1++) {
      part1 = &p_part1[p1];
      printf("%d %f %f %f %d %d %d\n", part1->r.identity,
	     part1->r.p[0], part1->r.p[1], part1->r.p[2],
	     part1->i[0], part1->i[1], part1->i[2]
	     );
    }
  }
  */

  free(send_coords);
  free(send_forces);
  free(recv_coords);
  free(recv_forces);
  free(sizes);
}

#endif
