#include <math.h>
#include "config.h"
#include "mmm2d.h"
#include "mmm-common.h"
#include "utils.h"
#include "specfunc.h"

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

MMM2D_struct mmm2d_params = { };

/****************************************
 * LOCAL FUNCTIONS
 ****************************************/

/* Bessel evaluation */
static void prepareBesselCutoffs(int P);

/* complex evaluation */
static void prepareBernoulliNumbers(int nmax);

/****************************************
 * NEAR FORMULA
 ****************************************/

/* cutoff error setup */
static void MMM2DSetPairwiseError(double error)
{
  int P, n, i;
  double T, err;
  double rho2m2max;

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
  T = log(part_error/16/M_SQRT2);
  for (i = 0; i < COMPLEX_STEP; i++)
    complexCutoff[i] = (int)ceil(T/log((i+1)/COMPLEX_FAC));
  prepareBernoulliNumbers(complexCutoff[COMPLEX_STEP-1]);

  /* polygamma, determine order */
  n = 0;
  rho2m2max = 1.0;
  do {
    create_mod_psi_up_to(n+1);

    err = fabs(2*mod_psi_even(n,rho2m2max));
    rho2m2max *= 0.5;
    n++;
  }
  while (err > 0.1*mmm2d_params.maxPWerror);
}

static void prepareBesselCutoff(int P)
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
  double pref;

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
  
  pref = C_2PISQR;
  for(l = 1; (l <= bon_order) && (l < 15); l++) {
    bon.e[l-1] = pref*bon_table[l]/l;
    pref *= C_2PISQR;
  }
  for (; l <= bon_order; l++) {
    if (l & 1)
      bon.e[l-1] =  2.0/l;
    else
      bon.e[l-1] = -2.0/l;      
  }
}

static void MMM2DcalcFromCopies(double *_pt, double *_Fx, double *_Fy,
				double *_Fz, double x, double y, double z)
{
  double rho_s = y*y + z*z;
  double pt, Fx, Fy, Fz;

  /* explicitly added potentials r_{-1,0} and r_{1,0} and constant
     part of potential */
  {
    register double cx    = x + 1;
    register double rinv2 = 1.0/(cx*cx + rho_s);
    register double rinv  = sqrt(rinv2);
    register double rinv3 = rinv*rinv2;
    pt = C_2LOG4PI + rinv;
    Fx = cx*rinv3;
    Fy =  y*rinv3;
    Fz =  z*rinv3;

    cx   = x - 1;
    rinv2 = 1.0/(cx*cx + rho_s);
    rinv  = sqrt(rinv2);
    rinv3 = rinv*rinv2;
    pt += rinv;
    Fx += cx*rinv3;
    Fy +=  y*rinv3;
    Fz +=  z*rinv3;
  }

  /* Bessel sum */
  /* the code is a bit lengthy since we two different chebychev expansions
     for the Bessel functions are used, one for 2<=arg<=8 and one for the
     rest (arg >= 2 pi)
  */
  {
    register int p, l;
    for (p = 1; p < besselEnd; p++) {
      register double k0Sum = 0;
      register double k1ySum = 0, k1Sum = 0;
      register double freq;
      register double rho, arg, ypl;
      double k0, k1;

      freq = C_2PI*p;

      /* series for 2 <= arg <= 8 */
      for (l = 1; l < besselCutoff[p-1]; l++) {
	ypl = y + l;
	rho = sqrt(ypl*ypl + z*z);
	arg = freq*rho;
	if (arg > 8) {
	  K0_1_i(arg, &k0, &k1);
	  k1 /= rho;
	  k0Sum  += k0;
	  k1ySum += k1*ypl;
	  k1Sum  += k1;
	  break;
	}
	K0_1_8(arg, &k0, &k1);
	k1 /= rho;
	k0Sum  += k0;
	k1ySum += k1*ypl;
	k1Sum  += k1;
      }
      /* rest with arg > 8 */
      for (l++; l < besselCutoff[p-1]; l++) {
	ypl = y + l;
	rho = sqrt(ypl*ypl + z*z);
	arg = freq*rho;
	K0_1_i(arg, &k0, &k1);
	k1 /= rho;
	k0Sum  += k0;
	k1ySum += k1*ypl;
	k1Sum  += k1;
      }

      /* series for 2 <= arg <= 8, negative l.
	 These are two independent loops since the switch occurs
	 at different l's.
      */
      for (l = 1; l < besselCutoff[p-1]; l++) {
	ypl = y - l;
	rho = sqrt(ypl*ypl + z*z);
	arg = freq*rho;
	if (arg > 8) {
	  K0_1_i(arg, &k0, &k1);
	  k1 /= rho;
	  k0Sum  += k0;
	  k1ySum += k1*ypl;
	  k1Sum  += k1;
	  break;
	}
	K0_1_8(arg, &k0, &k1);
	k1 /= rho;
	k0Sum  += k0;
	k1ySum += k1*ypl;
	k1Sum  += k1;
      }

      /* rest with arg > 8 */
      for (l++; l < besselCutoff[p-1]; l++) {
	ypl = y - l;
	rho = sqrt(ypl*ypl + z*z);
	arg = freq*rho;
	K0_1_i(arg, &k0, &k1);
	k1 /= rho;
	k0Sum  += k0;
	k1ySum += k1*ypl;
	k1Sum  += k1;
      }

      k0 = cos(freq*x);
      k1 = sin(freq*x);
      pt +=        k0*k0Sum;
      Fx +=   freq*k1*k0Sum;
      Fy +=   freq*k0*k1ySum;
      Fz += z*freq*k0*k1Sum;
    }
  }

  /* complex sum */
  {
    register double zeta_r, zeta_i;
    register double zet2_r, zet2_i;
    register double ztn_r,  ztn_i;
    register int end, n;

    ztn_r = zeta_r = z;
    ztn_i = zeta_i = y;
    zet2_r = zeta_r*zeta_r - zeta_i*zeta_i;
    zet2_i = 2*zeta_r*zeta_i;

    end = complexCutoff[(int)ceil(COMPLEX_FAC*rho_s)];
    for (n = 0; n < end; n++) {
      register double tmp_r;
      Fy -= 2*(n + 1)*bon[n]*ztn_i;
      Fz += 2*(n + 1)*bon[n]*ztn_r;
      pt -= bon[n]*(ztn_r*zeta_r - ztn_i*zeta_i);

      tmp_r = ztn_r*zet2_r - ztn_i*zet2_i;
      ztn_i = ztn_r*zet2_i + ztn_i*zet2_r;
      ztn_r = tmp_r;
    }
  }

  /* psi sum */
  {
    register int n;
    register double rho_snm1; /* rho^{2n-2} */

    /* n = 0 inflicts only Fx and pot */
    pt -= mod_psi_even(0, x);
    Fx += mod_psi_odd(0, x);

    rho_snm1 = 1.0;
    for (n = 1;; n++) {
      register double deriv  = 2*n;
      register double mpe    = mod_psi_even(n, x);
      register double mpo    = mod_psi_odd(n, x);
      register double rho_sn = rho_snm1*rho_s;

      pt -=       rho_sn  *mpe;
      Fx +=       rho_sn  *mpo;
      Fy += deriv*rho_snm1*mpe*y;
      Fz += deriv*rho_snm1*mpe*z;

      if (fabs(deriv*rho_snm1*mpe) < part_error)
	break;

      rho_snm1 = rho_sn;
    }
  }
  *_pt = pt;
  *_Fx = Fx;
  *_Fy = Fy;
  *_Fz = Fz;
}

double MMM2DselfEnergy()
{
  double pt = 0;
  double dummy;
  MMM2DcalcFromCopies(&pt, &dummy, &dummy, &dummy, 0, 0, 0);
  return pt;
}

void add_mmm2d_pair_force(Particle *p1, Particle *p2, int type_num)
{
  register double rinv2;
  register double rinv;
  register double rinv3;

  x -= rint(x);
  y -= rint(y);

  MMM2DcalcFromCopies(pt, Fx, Fy, Fz, x, y, z);

  rinv2 = 1.0/(x*x + y*y + z*z);
  rinv  = sqrt(rinv2);
  rinv3 = rinv*rinv2;
  *pt += rinv;
  *Fx += rinv3*x;
  *Fy += rinv3*y;
  *Fz += rinv3*z;
}

void MMM2DdoInternalChecks()
{
  FILE *out;
  int i, j, n;
  double Wn, jfac;

  out = fopen("bon.log", "w");
  fprintf(stderr, "generating bon.log for %d bernoulli numbers\n",
	  bonMax);
  Wn = C_2PI*C_2PI;
  jfac = 2;
  for (j = 1; j <= bonMax; j++) {
    fprintf(out, "%d %f %f\n", j, bon[j-1], jfac*bon[j-1]/(2*Wn)*(2*j));
    Wn *= C_2PI*C_2PI;
    jfac *= (2*j+1)*(2*j+2);
  }
  fclose(out);

  out = fopen("bessel.log", "w");
  fprintf(stderr, "generating bessel.log for a cutoff of %d\n", besselEnd);
  for (j = 1; j < besselEnd; j++)
    fprintf(out, "%d %d\n", j, besselCutoff[j - 1]);
  fclose(out);

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
}

#endif
