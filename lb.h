/*****************************************************************************/
/*                     basic types for machine independence                  */
/*****************************************************************************/

#ifndef LB_H
#define LB_H

typedef signed int         T_INT32;
typedef unsigned int       T_UINT32;
typedef unsigned int       T_UINT;

typedef float              T_SINGLE;
typedef signed char        T_BYTE;
typedef unsigned char      T_UBYTE;
typedef int                T_BOOL; 


typedef char              *TP_CHAR;
typedef T_BYTE            *TP_BYTE;
typedef T_UBYTE           *TP_UBYTE;
typedef int               *TP_INT; 
typedef T_UINT            *TP_UINT; 
typedef T_INT32           *TP_INT32; 
typedef T_UINT32          *TP_UINT32; 
typedef double            *TP_DOUBLE;
typedef T_SINGLE          *TP_SINGLE;


#ifndef FALSE
#define FALSE      ((T_BOOL)0)
#endif

#ifndef TRUE
#define TRUE       ((T_BOOL)-1)
#endif


#define MAXCOMTYPE    4              /* max type of commands and max # of commands */
#define MAXCOMMAND    20
#define SPACE_DIM     3
#define D_SPACE_DIM   3.0
#define SQRT2         1.4142135624
#define INVSQRT2      0.70710678119
#define SPACE_DIM_INV 0.33333333333333
#define MASSRATIO     0.5           /* mass ratio of fluid particles and monomers (fluid/mono) */
#define max_corr      900           /* max. correlation time for fluid correlations */ 
#define max_n_veloc   30            /* not more than 30 velocities for particles to take */
#define MAX_CHAIN     1200         /* max. length of one chain */
#define MAX_FOURIER   400           /* max. of k-vectors in StrFct calc. */
#define MAX_FOUR      40           /* max. of k-vectors in Fourier transform. */
#define N_MODE        3             /* number k-values for Fourier transform */
#define BOTT          1e-10         /* for comparison with DOUBLE zero */   
#define MAX_RANREPEATE 3


/* data structure definitions for whole program */

typedef int      T_IVECTOR [SPACE_DIM];
typedef double   T_DVECTOR [SPACE_DIM];
typedef double   T_DMATRIX [SPACE_DIM][SPACE_DIM];

typedef T_DVECTOR *TP_DVECTOR;
typedef T_IVECTOR *TP_IVECTOR;
typedef T_DMATRIX *TP_DMATRIX;

typedef T_INT32  T_NEIGHPE [2*SPACE_DIM];
typedef T_NEIGHPE *TP_NEIGHPE;

typedef T_INT32  T_GRIDPOS [SPACE_DIM+1];
typedef T_GRIDPOS *TP_GRIDPOS;

typedef double T_RANSTORE [MAX_RANREPEATE][SPACE_DIM];
typedef T_RANSTORE *TP_RANSTORE; 

typedef int T_CELLNEIGH [13];   
typedef T_CELLNEIGH *TP_CELLNEIGH;



/* UNITS:
 *  
 * I. lattice part:
 *
 * length scale is fixed by a which is the grid distance (in LJ units)           *
 * time   scale is fixed by tau which is the time step for propagation           *
 * mass   of particles is always one                                             *
 * energy scale is then given as a/tau^2                                         *
 *                                                                               *
 * in the lattice code all calculations are done with a and tau equal one,	 *
 * the end results are then multiplied by suitable a and tau powers		 *
 * to give the real world result and vice versa						 *
 * for convenience, after all physical variables, the factor for transforming    *
 * from lattice boltzmann to physical units is given.                            *

 * II. polymer part:
 *
 * All stuff is in LJ units which is therefore what we called the real world     *
 * above. in the coupling routine, we transform first the lattice results to     *
 * LJ units and then backwards when calculating the fluid particle shift         *
 * 
 * By this procedure each part can work in units where most used constants r ONE */

typedef struct {
  int     gridpoints;   /* # of gridpoints in each dimension  */
  double  rho;          /* number density [1/a^3] */
  double  viscos;      /* kinematic viscosity [a^2/tau] */
  double  friction;
  double  agrid;        /* grid distance [LJ units] */
  double  tau;         /* time step to propagate from one gridpoint to next [LJ] also used for Verlet step (preleminary */

  double  lambda;      /* relaxation time [tau] for non-eq distribution, not independent from viscos */
  double  c_sound_sq;  /* (sound velocity)^2 [a^2/tau^2]  */
  double  current_type;     /* current density [1/a^2/tau] */
  char    currentpar[200];
  T_INT32 ranseed;     /* seed for random number generator   */
  T_BOOL  boundary;    /* boundaries yes/no */
  int     chainlen;        
          
} LB_structure;

extern LB_structure compar; 

typedef struct {

  int offset;
  int doffset;
  int stride;
  int skip;
  int nblocks;
} T_PLANEINFO;

typedef T_PLANEINFO *TP_PLANEINFO;

/** If non-zero, the lb populations will no be updated during force recalc. */
extern int transfer_momentum;

void halo_init();

/*
** init halo stuff (plane information)
**/


void halo_update(); 
/* 
** update halo region
**/


void finishMPI();

typedef enum  {TERM_continue = 0,
	       TERM_kill     = 1, 
	       TERM_cpulimit = 2,
	       TERM_shutdown = 3,
	       TERM_sysstop  = 4,
	       TERM_error    = 5,
	       TERM_ok       = 6} T_TERM_FLAG;


void  InstallSignalHandler (void);
/*
** Install the signal handler for the program, which captures
** SIG_STOP, SIG_CPULIM, SIG_SHUTDN signal and checks for
** the existance of the 'sysStop' file.
** If one of these conditions are fullfilled the handler
** sets the terminate flag.
*/

T_TERM_FLAG Terminate (void);
/*
** Returns the status of the terminate flag
*/
		  
void FatalError (int numPe, char* where, char* what);
/*
** Prints an error message and exits the program
*/

void RecoverableError (int numPe, char* where, char* what);
/*
** Prints an error message
*/

void Warning (int numPe, char* where, char* what);
/*
** Prints a warning message
*/

void ExitProgram (T_TERM_FLAG termFlag);
/*
** Stop the program, showing a message according to the status of 'termFlag'
*/

void LB_Run();

/*************************** some simple macros ******************************/

#define sqr(a)         ((a)*(a))
#define cube(a)        ((a)*(a)*(a))
#define quad(a)        (((a)*(a))*((a)*(a)))

#ifndef max
#define max(a,b)       ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)       ((a) < (b) ? (a) : (b))
#endif

#define max3(a,b,c)    (max(max(a,b),c))
#define min3(a,b,c)    (min(min(a,b),c))

void calc_fluid_chain_interaction(int iamghost);
int lb_callback(Tcl_Interp *interp, void *_data);

void calc_lbforce();
void LB_propagate();
int lb_parser(Tcl_Interp * interp, int argc, char ** argv);
int printlbIAToResult(Tcl_Interp *interp);
void thermo_init_lb();

#endif /* LB_H */









