/*------------------------------------------------------------

  SUBUNIT:  p3m_v2.c

  NEEDS:    mymath.h, mymath.c, fixed_constants.h, p3m_v2.h
            ...and DXML!

  PURPOSE:  performs the particle-mesh calculations of the 
            P3M-method described in

            R.W. Hockney and J.W. Eastwood: Computer Simulation 
            Using Particles, IOP Bristol 1988

            see also: How to mesh up Ewald sums (I): A theoretical and
	    numerical comparison of various particle mesh routines,
	    M. Deserno, C. Holm, J. Chem. Phys. Vol 109, 7678 (1998)

  COMMENTS: - Contains a net charge correction in case the 
	      system is not electrostatically neutral.

  VERSION:  20 January 1999

  AUTHOR:   Markus Deserno

  CONTENT:  P3M_init:
              Initializes the P3M package by doing the following
	      things: Import the quantities systemlength, particle
	      number, alpha, Mesh, charge assignment order, Coulomb
	      prefactor. Initializes the
	      DXML-3d-FFT-routine. Interpolate the charge assignment
	      function. Calculate the shifted mesh indices. Calculate
	      the differential operator (-> ik-differentiation). And
	      last but not least: calculate the optimal influence
	      function of Hockney/Eastwood.
	      If any of the parameters important for this influence
	      function is changed during the simulation (e.g., the box
	      length), P3M must be exited and re-initialized with the
	      new parameters.  NOTE: Although the influence function
	      is independent of the number of particles, this number
	      is transported to the P3M unit only through the function
	      P3M_init. Thus, if this number changes, P3M has to be
	      reinitialized as well.
	    P3M_dipole:
	      Computes the dipole contribution to the electrostatic energy 
	      and forces in the Ewald method. This is done by using the
	      UNFOLDED particle coordinates, as is suggested in:
	      Jean-Michel Caillol: Comments on the numerical simulations of
	      electrolytes in periodic boundary conditions, J. Chem. Phys.,
	      Vol. 101, No. 7, 1. October 1994.
            P3M_perform:
              Performs one particle mesh calculation. As its input it
              needs (pointers to) the arrays of the particle
              coordinates, charges and forces. Also one can set two
              flags, which decide, whether or not the Coulomb
              forces/Energy should be calculated.
            P3M_perform_subset:
	      Calculates the k-space contribution to the energy coming from
	      interactions of particles within the subset {start,...,end}.
            P3M_exit:
              Exits P3M.

  USE:      The four functions P3M_init, P3M_exit and P3M_perform and
            P3M_perform_subset have the following calling syntax:
	    void P3M_init(double [length of the system], 
                          int    [number of particles],
                          p3m_struct *[pointer on p3m_struct, defined in p3m_v2.h], 
	    void P3M_exit(void);
	    double P3M_dipole(double [][3] particle positions, 
			      double *[array of particle charges],		  
			      double [][3] particle forces,
                             Returns the dipole energy.
	    void P3M_perform(double  [][3] particle positions,
			     double  *[array of particle charges], 
		             double  [][3] particle forces,
			     double  *[(k-space contribution) Coulomb Energy]);
	    double P3M_perform_subset(double  [][3] particle positions,
				      double  *[array of particle charges], 
				      int     [first particle in subset], 
				      int     [last particle in subset]);
				     Returns the energy.
				      
  CHANGES:  - (19.11.98) coordinates and forces are now matrices
            - (19.01.99) phi_{x,y,z}_{re,im} has been redefined
            - (20.01.99) Dn[mesh/2] has been set to 0 for symmetry reasons

------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

#include <dxmldef.h> 
#include "global.h"

#include "utils.h"
#include "global.h"
#include "particle_data.h"
/*#include <fixed_constants.h>*/
#include "p3m_parallel.h"

int zfft_apply_3d_( char *, char *, char *, ... );
int zfft_init_3d_( int *, int *, int *, DXML_Z_FFT_STRUCTURE_3D *, int * );
int zfft_exit_3d_( DXML_Z_FFT_STRUCTURE_3D * ); 
/*----------------------------------------------------------------------*/

/* some variables, internal to p3m_v2.c, which are shared among the
   routines needed for P3M */

static double  L,L2,L3,Li; /* box length, its square, cube and inverse */
static int     NP;         /* number of particles */
static double  alpha;      /* Ewald parameter */
static int     mesh;       /* number of mesh points */
static double  dmesh;      /* dito, but casted to double */
static int     P;          /* charge assignment order */
static double  prefactor;  /* prefactor of the Coulomb potential */
static double  surreps;    /* epsilon off the surroundings during L -> infinity */
p3m_struct p3m;

/* particle-coordinates, folded back into the box */
double ** cooP;

/* stores the meshpoints, shifted by mesh/2 */
static double  meshift[MAXMESH];

/* global index of a charged particle */
int     *global;

/* reference points in the lattice, 
   needed for the charge assignment */
int     **G;

/* stores the Fourier transformed differential 
   operator, however, this is done on the level 
   of n-vectors and not k-vectors. */
static double  Dn[MAXMESH];

/* charge of the charged particles */
double  *QP;

/* mesh contributions of the charged particles */
double  *QL;

/* stores the influence function: */
static double  Ghat[MAXMESH][MAXMESH][MAXMESH];
/* stores the interpolation of the charge assignment functions */
static double  intCAF[P_max][2*MaxInterpol+1];

/* ...and for the FFT: */
static double  Q_re[MAXMESH][MAXMESH][MAXMESH];
static double  Q_im[MAXMESH][MAXMESH][MAXMESH];
static double  phi_x_re[MAXMESH][MAXMESH][MAXMESH];
static double  phi_x_im[MAXMESH][MAXMESH][MAXMESH];
static double  phi_y_re[MAXMESH][MAXMESH][MAXMESH];
static double  phi_y_im[MAXMESH][MAXMESH][MAXMESH];
static double  phi_z_re[MAXMESH][MAXMESH][MAXMESH];
static double  phi_z_im[MAXMESH][MAXMESH][MAXMESH];

/* init-structure for the FFT: */
DXML_Z_FFT_STRUCTURE_3D fft_struct;
/* some funny FFT parameter, see DXML-manual */
static int     fft_ni_stride = 1;

/*----------------------------------------------------------------------*/

/* internal functions */

static void    perform_aliasing_sums(int NX, int NY, int NZ, 
				     double *nominatorX, double *nominatorY, 
				     double *nominatorZ, double *denominator);
static void    calculate_differential_operator(void);
static void    calculate_influence_function(void);
static void    interpolate_charge_assignment_function(void);
static void    calculate_meshift(void);
static double  sinc(double d);

/*----------------------------------------------------------------------*/

/* external functions */

void   P3M_init(double *length, int particlenumber, p3m_struct *p3m);
void   P3M_exit(void);
double P3M_dipole(int force);
void   P3M_perform(int force, double *E_Coulomb_P3M);
double P3M_perform_subset(int start, int end);
     
/*----------------------------------------------------------------------*/

static void perform_aliasing_sums(int NX, int NY, int NZ, 
				  double *nominatorX, double *nominatorY, 
				  double *nominatorZ, double *denominator)
{
  /* calculates the aliasing sums in the nominator and denominator of the 
     expression for the optimal influence function (see  Hockney/Eastwood: 
     8-22, p. 275).
     NX, NY, NZ: components of the n-vector for which the aliasing sum is 
                 to be performed.
     *nominatorX,*nominatorY,*nominatorZ : x-, y-, and z-component of the 
                                           aliasing sum.
     *denominator : aliasing sum in the denominator. */
  
  double S1,S2,S3;
  double fak1,fak2,fak3;
  int    MX,MY,MZ;
  double NMX,NMY,NMZ;
  double NM2;
  double expo;

  static double exponent_limit = 30;

  fak1 = 1.0/dmesh;
  fak2 = SQR(PI/(alpha*L));

  *nominatorX = *nominatorY = *nominatorZ = *denominator = 0.0;

  for (MX = -Brillouin; MX <= Brillouin; MX++) {
    NMX = meshift[NX] + dmesh*MX;
    S1  = pow(sinc(fak1*NMX),2.0*P);
    for (MY = -Brillouin; MY <= Brillouin; MY++) {
      NMY = meshift[NY] + dmesh*MY;
      S2  = S1*pow(sinc(fak1*NMY),2.0*P);
      for (MZ = -Brillouin; MZ <= Brillouin; MZ++) {
	NMZ = meshift[NZ] + dmesh*MZ;
	S3  = S2*pow(sinc(fak1*NMZ),2.0*P);
	
	*denominator += S3;

	NM2 = SQR(NMX) + SQR(NMY) + SQR(NMZ);
	
	fak3 = ((expo=fak2*NM2)<exponent_limit) ? S3*exp(-expo)/NM2 : 0.0;

	*nominatorX += fak3 * NMX;
	*nominatorY += fak3 * NMY;
	*nominatorZ += fak3 * NMZ;
      }
    }
  }
}

/*----------------------------------------------------------------------*/

static void calculate_differential_operator(void)
{
  /* Calculates the Fourier transformed differential operator, 
     however, this is done on the level of n-vectors and not 
     k-vectors, i.e. the prefactor i*2*PI/L is missing! */
  
  int i;
  
  fprintf(stderr," - calculating differential operator");

  for (i=0; i<mesh; i++) 
    Dn[i] = (double)i - dround((double)i/dmesh)*dmesh;

  /* for symmetry reasons it is clever to set 
     Dn[mesh/2] to 0, for otherwise the expected 
     symmetry Dn[n] = -Dn[mesh-n] does not hold. */
  Dn[mesh/2] = 0.0;

  fprintf(stderr,"\n");
}

/*----------------------------------------------------------------------*/

static void calculate_influence_function(void)
{
  /* calculates the optimal influence function of Hockney and
     Eastwood, see: Hockney/Eastwood 8-22 (p275). Note the somewhat
     different convention for the prefactors, which is described in
     Deserno/Holm. */

  int    NX,NY,NZ;
  double Dnx,Dny,Dnz,Dn2;
  double fak1,fak2,fak3;
  double nominatorX,nominatorY,nominatorZ,denominator;

  fprintf(stderr," - calculating influence function with parameters\n");
  fprintf(stderr,"   mesh=%d, P=%d, alpha=%lf, L=%lf, Brillouin=%d",
	  mesh,P,alpha,L,Brillouin);
  
  fak1  = dmesh*dmesh*dmesh * 2.0 / L2;
  fak2  = SQR(PI/(alpha*L));

  for (NX = 0; NX < mesh; NX++) 
    for (NY = 0; NY < mesh; NY++) 
      for (NZ = 0; NZ < mesh; NZ++) 
	if ( (NX == 0) && (NY == 0) && (NZ == 0) )
	  Ghat[NX][NY][NZ]=0.0; 
	else {
	  perform_aliasing_sums(NX,NY,NZ,&nominatorX,&nominatorY,&nominatorZ,&denominator);

	  Dnx = Dn[NX]; 
	  Dny = Dn[NY]; 
	  Dnz = Dn[NZ]; 

	  Dn2 = SQR(Dnx) + SQR(Dny) + SQR(Dnz);
	  
	  if ( Dn2 > 1e-7 ) {
	    fak3  = Dnx*nominatorX + Dny*nominatorY + Dnz*nominatorZ; 
	    fak3 /= ( Dn2 * SQR(denominator) ); 
	    Ghat[NX][NY][NZ] = fak3*fak1; 
	  } 
	  else Ghat[NX][NY][NZ] = 0.0;
	}
  
  fprintf(stderr,"\n");
}

/*----------------------------------------------------------------------*/

static void interpolate_charge_assignment_function(void)
{
  /* Interpolates the P-th order charge assignment function from
     Hockney/Eastwood 5-189 (or 8-61). The following charge fractions
     are also tabulated in Deserno/Holm. */

  double dInterpol=(double)MaxInterpol, x;
  long   i;

  fprintf(stderr," - interpolating the order-%d charge assignment function",P);

  switch (P) {
  case 1 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 1.0;
    }
  } break;
  case 2 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 0.5-x;
      intCAF[1][i+MaxInterpol] = 0.5+x;
    }
  } break;
  case 3 : { 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = 0.5*SQR(0.5 - x);
      intCAF[1][i+MaxInterpol] = 0.75 - SQR(x);
      intCAF[2][i+MaxInterpol] = 0.5*SQR(0.5 + x);
    }
  } break;
  case 4 :{ 
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = ( 1.0+x*( -6.0+x*( 12.0-x* 8.0)))/48.0;
      intCAF[1][i+MaxInterpol] = (23.0+x*(-30.0+x*(-12.0+x*24.0)))/48.0;
      intCAF[2][i+MaxInterpol] = (23.0+x*( 30.0+x*(-12.0-x*24.0)))/48.0;
      intCAF[3][i+MaxInterpol] = ( 1.0+x*(  6.0+x*( 12.0+x* 8.0)))/48.0;
    }
  } break;
  case 5 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (  1.0+x*( -8.0+x*(  24.0+x*(-32.0+x*16.0))))/384.0;
      intCAF[1][i+MaxInterpol] = ( 19.0+x*(-44.0+x*(  24.0+x*( 16.0-x*16.0))))/ 96.0;
      intCAF[2][i+MaxInterpol] = (115.0+x*       x*(-120.0+x*       x*48.0))  /192.0;
      intCAF[3][i+MaxInterpol] = ( 19.0+x*( 44.0+x*(  24.0+x*(-16.0-x*16.0))))/ 96.0;
      intCAF[4][i+MaxInterpol] = (  1.0+x*(  8.0+x*(  24.0+x*( 32.0+x*16.0))))/384.0;
    }
  } break;
  case 6 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (  1.0+x*( -10.0+x*(  40.0+x*( -80.0+x*(  80.0-x* 32.0)))))/3840.0;
      intCAF[1][i+MaxInterpol] = (237.0+x*(-750.0+x*( 840.0+x*(-240.0+x*(-240.0+x*160.0)))))/3840.0;
      intCAF[2][i+MaxInterpol] = (841.0+x*(-770.0+x*(-440.0+x*( 560.0+x*(  80.0-x*160.0)))))/1920.0;
      intCAF[3][i+MaxInterpol] = (841.0+x*(+770.0+x*(-440.0+x*(-560.0+x*(  80.0+x*160.0)))))/1920.0;
      intCAF[4][i+MaxInterpol] = (237.0+x*( 750.0+x*( 840.0+x*( 240.0+x*(-240.0-x*160.0)))))/3840.0;
      intCAF[5][i+MaxInterpol] = (  1.0+x*(  10.0+x*(  40.0+x*(  80.0+x*(  80.0+x* 32.0)))))/3840.0;
    }
  } break;
  case 7 : {
    for (i=-MaxInterpol; i<=MaxInterpol; i++) {
      x=i/(2.0*dInterpol);
      intCAF[0][i+MaxInterpol] = (    1.0+x*(   -12.0+x*(   60.0+x*( -160.0+x*(  240.0+x*(-192.0+x* 64.0))))))/46080.0;
      intCAF[1][i+MaxInterpol] = (  361.0+x*( -1416.0+x*( 2220.0+x*(-1600.0+x*(  240.0+x*( 384.0-x*192.0))))))/23040.0;
      intCAF[2][i+MaxInterpol] = (10543.0+x*(-17340.0+x*( 4740.0+x*( 6880.0+x*(-4080.0+x*(-960.0+x*960.0))))))/46080.0;
      intCAF[3][i+MaxInterpol] = ( 5887.0+x*          x*(-4620.0+x*         x*( 1680.0-x*        x*320.0)))   /11520.0;
      intCAF[4][i+MaxInterpol] = (10543.0+x*( 17340.0+x*( 4740.0+x*(-6880.0+x*(-4080.0+x*( 960.0+x*960.0))))))/46080.0;
      intCAF[5][i+MaxInterpol] = (  361.0+x*(  1416.0+x*( 2220.0+x*( 1600.0+x*(  240.0+x*(-384.0-x*192.0))))))/23040.0;
      intCAF[6][i+MaxInterpol] = (    1.0+x*(    12.0+x*(   60.0+x*(  160.0+x*(  240.0+x*( 192.0+x* 64.0))))))/46080.0;
    }
  } break;
  default :{
    fprintf(stderr,"Error in function 'interpolate_charge_assignment_function':");
    fprintf(stderr,"Charge assignment order %d unknown.\nProgram terminated.\n\n",P);
    exit(1);
  }
  }

  fprintf(stderr,"\n");
}

/*----------------------------------------------------------------------*/

static void calculate_meshift(void)
{
  /* shifts the mesh points by mesh/2 */
  
  int i;

  fprintf(stderr," - calculating mesh-shift");
  
  for (i=0; i<mesh; i++) meshift[i] = i - dround(i/dmesh)*dmesh; 
  
  fprintf(stderr,"\n");
}

/*----------------------------------------------------------------------*/

double P3M_dipole(int force)
{
  /* Computes the dipole contribution to the electrostatic energy 
     and forces in the Ewald method. This is done by using the
     UNFOLDED particle coordinates, as is suggested in:
     Jean-Michel Caillol: Comments on the numerical simulations of
     electrolytes in periodic boundary conditions, J. Chem. Phys.,
     Vol. 101, No. 7, 1. October 1994.

     coo[][]    particle coordinates
     Q[]        particle charges
     For[][]    forces on the particles. If the user 
                supplies 0, the forces are not computed. */

  int    i;
  double sum[3] = { 0.0, 0.0, 0.0 };
  double E;
  double zero = 1e-5;
  double pref = prefactor * 4.0 * PI / ((1.0 + 2.0*surreps) * L3);

  /* dipole moment: */
  for (i = 0; i < NP; i++) if (fabs(particles[i].q) > zero) {
    sum[0] += particles[i].q * particles[i].p[0];
    sum[1] += particles[i].q * particles[i].p[1];
    sum[2] += particles[i].q * particles[i].p[2];
  }

  E = 0.5 * pref * ( SQR(sum[0]) + SQR(sum[1]) + SQR(sum[2]) );
  
  sum[0] *= pref;
  sum[1] *= pref;
  sum[2] *= pref;
  
  if (force) 
    for (i = 0; i < NP; i++) 
      if (fabs(particles[i].q) > zero) {
	particles[i].f[0] -= particles[i].q * sum[0];
	particles[i].f[1] -= particles[i].q * sum[1];
	particles[i].f[2] -= particles[i].q * sum[2];
      }

  return E;
}

/*----------------------------------------------------------------------*/

void P3M_perform(int force, double *E_Coulomb_P3M)
{
  /* Calculates the particle-mesh-contribution of the P3M routine.

     coo[][]       : coordinates of the particles.
     Q[]           : charges of the particles (double).
     For[][]       : force-components of the particles. (0, if not desired)
     E_Coulomb_P3M : pointer to the k-space energy contribution of P3M. (0, if not desired)
  */
  
  /* counting variables: */
  int i, j, k, l, m; 
  /* variables needed for the FFT: */
  int Lda, Ldb, sx, sy, sz, status;   
  /* fast modulo, works only for mesh=2^...: */
  int MESHMASK;
  /* supporting variables */
  double d1,d2,Hi,MI2,modadd1,modadd2;
  /* accellerates the charge assignment */
  double T1,T2,T3;
  /* alternative reference to the arrays G[i][3] */
  int Gi0, Gi1, Gi2;
  /* arguments of the array intCAF */
  int xarg,yarg,zarg;
  /* positions on the lattice */
  int xpos,ypos,zpos;
  /* That much left to the reference points the charge
     assignment starts. There is also an additional
     summand 'mesh' for the modulo: */
  int assignshift;
  /* number of charged particles */
  int QZahl;
  /* SQUARE of SUM of charges */
  static double  sum_q_2; 
  /* SUM of SQUARE of charges */
  static double  sum_q2; 

  cooP      = malloc(sizeof(int *)*n_particles);
  for(i=0;i<n_particles;i++)
    cooP[i] = malloc(sizeof(int)*3);
  global = malloc(sizeof(int)*n_total_particles);
  G      = malloc(sizeof(int *)*n_particles);
  for(i=0;i<n_particles;i++)
    G[i] = malloc(sizeof(int)*3);
  QP = malloc(sizeof(double)*n_total_particles);
  QL = malloc(sizeof(double)*n_total_particles*(P_max*P_max*P_max));

  Lda = Ldb = MAXMESH; 
  sx = sy = sz = 1; 
  MESHMASK = mesh-1;     /* THIS REQUIRES MESH TO BE A POWER OF 2 !!! */
  Hi = (dmesh = (double)mesh) / L;
  MI2 = 2.0*(double)MaxInterpol;
  assignshift = mesh-(P-1)/2;  /* INTEGER division is vital here! */

  /* picking up the charged particles, folding them into the box and calcu-
     lating the square of the sum and the sum of the square of the charges: */
  QZahl = 0;
  sum_q_2 = sum_q2 = 0.0;
  for (i=0; i<n_particles; i++) if (fabs(particles[i].q)>1e-5) {
    cooP[QZahl][0] = particles[i].p[0] - floor(particles[i].p[0]*Li)*L;
    cooP[QZahl][1] = particles[i].p[1] - floor(particles[i].p[1]*Li)*L;
    cooP[QZahl][2] = particles[i].p[2] - floor(particles[i].p[2]*Li)*L;

    QP[QZahl] = particles[i].q;
    sum_q_2  += QP[QZahl];
    sum_q2   += SQR( QP[QZahl] );

    global[QZahl++] = i;
  }
  sum_q_2 *= sum_q_2;

  /********** CHARGE ASSIGNMENT **********/
  
  /* init Q_re and Q_im */
  for (i=0; i<mesh; i++) 
    for (j=0; j<mesh; j++) 
      for (k=0; k<mesh; k++)
	Q_re[i][j][k] = Q_im[i][j][k] = 0.0;

  /* prepare to distinguish between odd and even interpolation 
     order (note that modadd1+modadd2 >= 0): */
  switch (P) {
  case 2 : case 4 : case 6 : 
    { modadd1 = 0.5; modadd2 = -0.5;} break;
  case 1 :case 3 : case 5 : case 7 : 
    { modadd1 = 0.0; modadd2 =  0.5;} break;
  default : {
    fprintf(stderr,"Error in function 'P3M_perform':\n");
    fprintf(stderr,"Charge assignment order P=%d unknown.\nProgram terminated\n\n",P);
    exit(1);
  } break;
  }

  /* The following part looks slightly different depending 
     on whether or not one needs the forces:  */

  if (force) {   /* if forces should be computed... */
    /* Determine the reference points in the lattice as well 
       as the minimum distance to the nearest lattice point: */

    m = 0; /* <-- enumerate the charge fractions consecutively */
    for (i=0; i<QZahl; i++) {
      d1      = cooP[i][0]*Hi + modadd1; 
      G[i][0] = Gi0 = (int)(d1 + modadd2) + assignshift;
      xarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      
      d1      = cooP[i][1]*Hi + modadd1; 
      G[i][1] = Gi1 = (int)(d1 + modadd2) + assignshift;
      yarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      
      d1      = cooP[i][2]*Hi + modadd1; 
      G[i][2] = Gi2 = (int)(d1 + modadd2) + assignshift;
      zarg    = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      /* Note: The (int) in (int)(d1 + modadd2) is OK, since 
	 d1 + modadd2 = xP[i]*Hi + modadd1 + modadd2 >= 0,
	 thus one does not need the 'floor' function */
      
      /* Calculate the mesh based charges in Q_re[][][] 
	 and the charge fractions of each charge QL[]: */
      for (j = 0; j < P; j++) {
	xpos = (Gi0 + j) & MESHMASK;
	T1 = QP[i] * intCAF[j][xarg];
	for (k = 0; k < P; k++) {
	  ypos = (Gi1 + k) & MESHMASK;
	  T2 = T1 * intCAF[k][yarg];
	  for (l = 0; l < P; l++) {
	    zpos = (Gi2 + l) & MESHMASK;
	    T3 = T2 * intCAF[l][zarg];
	    
	    Q_re[xpos][ypos][zpos] += (QL[m++] = T3);
	  }
	}
      }
    }
  } else {   /* ...else if (For)... */
    for (i=0; i<QZahl; i++) {
      d1   = cooP[i][0]*Hi + modadd1; 
      Gi0  = (int)(d1 + modadd2) + assignshift;
      xarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      
      d1   = cooP[i][1]*Hi + modadd1; 
      Gi1  = (int)(d1 + modadd2) + assignshift;
      yarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      
      d1   = cooP[i][2]*Hi + modadd1; 
      Gi2  = (int)(d1 + modadd2) + assignshift;
      zarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
      
      /* Calculate the mesh based charges in Q_re[][][] 
	 and the charge fractions of each charge QL[]: */
      for (j = 0; j < P; j++) {
	xpos = (Gi0 + j) & MESHMASK;
	T1 = QP[i] * intCAF[j][xarg];
	for (k = 0; k < P; k++) {
	  ypos = (Gi1 + k) & MESHMASK;
	  T2 = T1 * intCAF[k][yarg];
	  for (l = 0; l < P; l++) {
	    zpos = (Gi2 + l) & MESHMASK;
	    T3 = T2 * intCAF[l][zarg];
	    
	    Q_re[xpos][ypos][zpos] += T3;
	  }
	}
      }
    }
  } /* ...end if (For)... */
  
  /* perform the forward FFT, using a DXML-routine: */
  status = zfft_apply_3d_("R","R","F", Q_re, Q_im, Q_re, Q_im, 
			  &Lda, &Ldb, &fft_struct, &sx, &sy, &sz); 

  /* If desired, calculate the k-space contribution to the Coulomb Energy */
  if (E_Coulomb_P3M) {    /* if energy should be computed... */
    *E_Coulomb_P3M = 0.0; /* Note that *E_Coulomb_P3M is first set to 0 */
    for (i=0; i<mesh; i++)
      for (j=0; j<mesh; j++)
	for (k=0; k<mesh; k++)
	  *E_Coulomb_P3M += Ghat[i][j][k] * ( SQR(Q_re[i][j][k]) + SQR(Q_im[i][j][k]) );
    *E_Coulomb_P3M *= ( prefactor * L / (dmesh*dmesh*dmesh*4.0*PI) );
    /* self energy and net charge correction: */
    *E_Coulomb_P3M -= prefactor * ( sum_q2 * alpha / wupi  +  
				    sum_q_2 * PI / (2.0*L3*SQR(alpha)) );
  }

  /* Only if the forces on the particles are needed, the following has to be done: */
  if (force) {   /* if forces should be computed... */
    /* Calculate the supporting arrays phi_?_??: */
    for (i=0; i<mesh; i++)
      for (j=0; j<mesh; j++)
	for (k=0; k<mesh; k++) {  
	  d1 = Ghat[i][j][k] * Q_re[i][j][k];
	  d2 = Ghat[i][j][k] * Q_im[i][j][k];
	  
	  /* old definition: */
	  /*
	    phi_x_re[i][j][k] = d1 * Dn[i];
	    phi_x_im[i][j][k] = d2 * Dn[i];
	    phi_y_re[i][j][k] = d1 * Dn[j];
	    phi_y_im[i][j][k] = d2 * Dn[j];
	    phi_z_re[i][j][k] = d1 * Dn[k];
	    phi_z_im[i][j][k] = d2 * Dn[k];
	  */
	  /* new definition: */
	  phi_x_re[i][j][k] = -d2 * Dn[i];
	  phi_x_im[i][j][k] =  d1 * Dn[i];
	  phi_y_re[i][j][k] = -d2 * Dn[j];
	  phi_y_im[i][j][k] =  d1 * Dn[j];
	  phi_z_re[i][j][k] = -d2 * Dn[k];
	  phi_z_im[i][j][k] =  d1 * Dn[k];
	}
    
    /* perform the backward FFT, using a DXML-routine: */
    status = zfft_apply_3d_("R","R","B", phi_x_re, phi_x_im, phi_x_re, phi_x_im, 
			    &Lda, &Ldb, &fft_struct, &sx, &sy, &sz); 
    status = zfft_apply_3d_("R","R","B", phi_y_re, phi_y_im, phi_y_re, phi_y_im, 
			    &Lda, &Ldb, &fft_struct, &sx, &sy, &sz); 
    status = zfft_apply_3d_("R","R","B", phi_z_re, phi_z_im, phi_z_re, phi_z_im, 
			    &Lda, &Ldb, &fft_struct, &sx, &sy, &sz); 

    m = 0; /* <-- enumerate the charge fractions consecutively */
    for (i = 0; i < QZahl; i++) {
      for (j = 0; j < P; j++) {
	xpos = (G[i][0] + j) & MESHMASK;
	for (k = 0; k < P; k++) {
	  ypos = (G[i][1] + k) & MESHMASK;
	  for (l = 0; l < P; l++) {
	    zpos = (G[i][2] + l) & MESHMASK;
	    
	    d1 = prefactor * QL[m++];
	    
	    /* old relation between phi and force */
	    /*
	      For[global[i]][0] += d1 * phi_x_im[xpos][ypos][zpos];
	      For[global[i]][1] += d1 * phi_y_im[xpos][ypos][zpos];
	      For[global[i]][2] += d1 * phi_z_im[xpos][ypos][zpos];
	    */

	    /* new relation between phi and force */
	    particles[global[i]].f[0] -= d1 * phi_x_re[xpos][ypos][zpos];
	    particles[global[i]].f[1] -= d1 * phi_y_re[xpos][ypos][zpos];
	    particles[global[i]].f[2] -= d1 * phi_z_re[xpos][ypos][zpos];
	    
	  }
	}  
      }
    }
  }
  free(cooP);
  free(global);
  free(G);
  free(QP);
  free(QL);
}

/*----------------------------------------------------------------------*/

double P3M_perform_subset(int start, int end)
{
  /* Returns the particle-mesh-contribution of the P3M routine to the 
     energy, if one restricts to particles within the subset {start,...,end}.
     coo[][] : coordinates of the particles.
     Q[]     : charges of the particles (double).
     start   : first particle in the subset.
     end     : last particle in the subset.
  */
  
  /* counting variables: */
  int i, j, k, l; 
  /* variables needed for the FFT: */
  int Lda, Ldb, sx, sy, sz, status;   
  /* fast modulo, works only for mesh=2^...: */
  int MESHMASK;
  /* supporting variables */
  double d1,d2,Hi,MI2,modadd1,modadd2;
  /* accellerates the charge assignment */
  double T1,T2,T3;
  /* alternative reference to the arrays G[i][] */
  int Gi0,Gi1,Gi2;
  /* arguments of the array intCAF */
  int xarg,yarg,zarg;
  /* positions on the lattice */
  int xpos,ypos,zpos;
  /* That much left to the reference points the charge
     assignment starts. There is also an additional
     summand 'mesh' for the modulo: */
  int assignshift;
  /* number of charged particles */
  int QZahl;
  /* SQUARE of SUM of charges */
  static double  sum_q_2; 
  /* SUM of SQUARE of charges */
  static double  sum_q2; 
  /* returned energy */
  double energy = 0.0;

  G = malloc(sizeof(int)+n_particles);
  for(i=0;i<n_particles;i++)
    G[i] = malloc(sizeof(int)*3);
  cooP = malloc(sizeof(int)+n_particles);
  for(i=0;i<n_particles;i++)
    cooP[i] = malloc(sizeof(int)*3);
  QP = malloc(sizeof(double)*n_total_particles);
  QL = malloc(sizeof(double)*n_total_particles*(P_max*P_max*P_max));

  Lda = Ldb = MAXMESH; 
  sx = sy = sz = 1; 
  MESHMASK = mesh-1;  /* THIS REQUIRES MESH TO BE A POWER OF 2 !!! */ 
  Hi = (dmesh = (double)mesh) / L;
  MI2 = 2.0*(double)MaxInterpol;
  assignshift = mesh-(P-1)/2;  /* INTEGER division is vital here! */

  /* picking up the charged particles, folding them into the box and calcu-
     lating the square of the sum and the sum of the square of the charges: */
  QZahl = 0;
  sum_q_2 = sum_q2 = 0.0;
  for (i = start; i <= end; i++) if (fabs(particles[i].q) > 1e-5) {  /* Note the '<='! */
    cooP[QZahl][0] = particles[i].p[0] - floor(particles[i].p[0]*Li)*L;
    cooP[QZahl][1] = particles[i].p[1] - floor(particles[i].p[1]*Li)*L;
    cooP[QZahl][2] = particles[i].p[2] - floor(particles[i].p[2]*Li)*L;

    QP[QZahl] = particles[i].q;
    sum_q_2  += QP[QZahl];
    sum_q2   += SQR( QP[QZahl] );
    
    global[QZahl++] = i;
  }
  sum_q_2 *= sum_q_2;

  /********** CHARGE ASSIGNMENT **********/
  
  /* init Q_re and Q_im */
  for (i = 0; i < mesh; i++) 
    for (j = 0; j < mesh; j++) 
      for (k = 0; k < mesh; k++)
	Q_re[i][j][k] = Q_im[i][j][k] = 0.0;

  /* prepare to distinguish between odd and even interpolation 
     order (note that modadd1+modadd2 >= 0): */
  switch (P) {
  case 2 : case 4 : case 6 : 
    { modadd1 = 0.5; modadd2 = -0.5;} break;
  case 1 :case 3 : case 5 : case 7 : 
    { modadd1 = 0.0; modadd2 =  0.5;} break;
  default : {
    fprintf(stderr,"Error in function 'P3M_perform':\n");
    fprintf(stderr,"Charge assignment order P=%d unknown.\nProgram terminated\n\n",P);
    exit(1);
  } break;
  }

  for (i = 0; i < QZahl; i++) {
    d1   = cooP[i][0]*Hi + modadd1; 
    Gi0  = (int)(d1 + modadd2) + assignshift;
    xarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
    
    d1   = cooP[i][1]*Hi + modadd1; 
    Gi1  = (int)(d1 + modadd2) + assignshift;
    yarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
    
    d1   = cooP[i][2]*Hi + modadd1; 
    Gi2  = (int)(d1 + modadd2) + assignshift;
    zarg = (int)( (d1 - dround(d1) + 0.5)*MI2 );
    
    /* Calculate the mesh based charges in Q_re[][][] 
       and the charge fractions of each charge QL[]: */
    for (j = 0; j < P; j++) {
      xpos = (Gi0 + j) & MESHMASK;
      T1 = QP[i] * intCAF[j][xarg];
      for (k = 0; k < P; k++) {
	ypos = (Gi1 + k) & MESHMASK;
	T2 = T1 * intCAF[k][yarg];
	for (l = 0; l < P; l++) {
	  zpos = (Gi2 + l) & MESHMASK;
	  T3 = T2 * intCAF[l][zarg];
	  
	  Q_re[xpos][ypos][zpos] += T3;
	}
      }
    }
  }
  
  /* perform the forward FFT, using a DXML-routine: */
  status = zfft_apply_3d_("R","R","F", Q_re, Q_im, Q_re, Q_im, 
			  &Lda, &Ldb, &fft_struct, &sx, &sy, &sz); 

  /* Calculate the k-space contribution to the Coulomb Energy */
  d1 = prefactor * L / (dmesh*dmesh*dmesh*4.0*PI);
  for (i = 0; i < mesh; i++)
    for (j = 0; j < mesh; j++)
      for (k = 0; k < mesh; k++) {
	d2 = ( SQR(Q_re[i][j][k]) + SQR(Q_im[i][j][k]) );
	energy += d1 * d2 * Ghat[i][j][k]; 
      }
  /* self energy and net charge correction: */
  energy -= prefactor * ( sum_q2 * alpha / wupi  +  
			  sum_q_2 * PI / (2.0*L3*SQR(alpha)) );
  
  free(cooP);
  free(global);
  free(G);
  free(QP);
  free(QL);
  return energy;
}

/*----------------------------------------------------------------------*/

void P3M_init(double *length, int particlenumber, p3m_struct *p3m)
{
  /* init the P3M-routine
     length           : length of the cubic system.
     particlenumber   : number of particles in the box.
     p3m              : pointer to p3m structure (defined in p3m_v2.h) */
  
  L         = length[0];
  L2        = L*L;
  L3        = L*L2;
  Li        = 1.0/L;
  NP        = particlenumber;
  alpha     = p3m->alpha;
  mesh      = p3m->mesh;
  dmesh     = (double)mesh;
  P         = p3m->P;
  prefactor = p3m->prefactor;
  surreps   = p3m->epsilon;

  if (mesh!=MAXMESH) {
    fprintf(stderr,"Desired mesh size (%d) and internal mesh "
	    "size (%d) are not equal!\n",mesh,MAXMESH);
    fprintf(stderr,"Program terminated!\n\n");
    exit(1);
  }
  
  /* init FFT: */
  zfft_init_3d_(&mesh,&mesh,&mesh,&fft_struct,&fft_ni_stride);

  interpolate_charge_assignment_function();
  calculate_meshift();
  calculate_differential_operator();
  calculate_influence_function();

  fprintf(stderr,"'P3M' successfully initialized!\n");
}

/*----------------------------------------------------------------------*/

void P3M_exit(void)
{
  /* deallocates the memory space needed for fft_struct */

  zfft_exit_3d_(&fft_struct);

  fprintf(stderr,"P3M deallocated!\n");
}

/*----------------------------------------------------------------------*/

static double sinc(double d)
{
  /* Calculates the sinc-function as sin(PI*x)/(PI*x) (same convention
     as in Hockney/Eastwood). In order to avoid divisions by 0,
     arguments, whose modulus is smaller than epsi, will be evaluated
     by an 8th order Taylor expansion of the sinc function. Note that
     the difference between sinc(x) and this expansion is smaller than
     0.235e-12, if x is smaller than 0.1. (The next term in the
     expansion is the 10th order contribution PI^10/39916800 * x^10 =
     0.2346...*x^12).  
     This expansion should also save time, since it reduces the number
     of function calls to sin().  */
  
  static double epsi =  0.1;

  static double   c2 = -0.1666666666667e-0;
  static double   c4 =  0.8333333333333e-2;
  static double   c6 = -0.1984126984127e-3;
  static double   c8 =  0.2755731922399e-5;

  double PId = PI*d, PId2;

  if (fabs(d) > epsi)
    return sin(PId) / PId;
  else {
    PId2 = SQR(PId);
    return 1.0 + PId2 * ( c2 + PId2 * ( c4 + PId2 * ( c6 + PId2 * c8 ) ) );
  }
  return 1.0;
}

/*----------------------------------------------------------------------*/
