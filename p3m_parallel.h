/*------------------------------------------------------------

HEADER-FILE fuer p3m_parallel.c

Letzte Aenderung: 22.02.2002

------------------------------------------------------------*/
/** \file p3m_parallel.h
    P3M header file for p3m_parallel.c
*/

#ifndef P3M_V2_H
#define P3M_V2_H

/* Maximal mesh size */
#define MAXMESH 32

/* Maximal charge assignment order */
#define P_max 7
                 
/* Maximal interpolation order */
#define MaxInterpol (4*50048)
   
/* number of Brillouin zones in the aliasing sums*/
#define Brillouin 1

/** \name p3m structure */
/*@}*/

/*----------------------------------------------------------*/

typedef struct {
  double alpha;     /** Ewald splitting parameter */
  int    mesh;      /** number of mesh points per coordinate direction */
  int    P;         /** charge assignment order */
  double rcut;      /** cutoff radius for real space electrostatics */
  double prefactor; /** Coulomb prefactor, e.g. Bjerrum*Temp */
  double epsilon;   /** epsilon of the "surrounding dielectric" */
}p3m_struct;

/*@}*/

/*----------------------------------------------------------*/

/** Initializes p3m */
/*@{*/
/*  Initializes the P3M package by doing the following
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
    reinitialized as well.*/
extern void   P3M_init(double *, int, p3m_struct *p3m);
/*@}*/

/** Exits p3m */
extern void   P3M_exit(void);

/** Computes dipole contribution */
/*@{*/
/*  Computes the dipole contribution to the electrostatic energy 
    and forces in the Ewald method. This is done by using the
    UNFOLDED particle coordinates, as is suggested in:
    Jean-Michel Caillol: Comments on the numerical simulations of
    electrolytes in periodic boundary conditions, J. Chem. Phys.,
    Vol. 101, No. 7, 1. October 1994.*/
extern double P3M_dipole(int force);
/*@}*/

/** Performs one particle mesh calculation*/
/*@{*/
/*  Performs one particle mesh calculation. As its input it
    needs (pointers to) the arrays of the particle
    coordinates, charges and forces. Also one can set two
    flags, which decide, whether or not the Coulomb
    forces/Energy should be calculated.*/
extern void   P3M_perform(int force, double *E_Coulomb_P3M);
/*@}*/

/** P3M_perform for only a subset */
/*@{*/
/* Calculates the k-space contribution to the energy coming from
   interactions of particles within the subset {start,...,end}.*/
extern double P3M_perform_subset(int, int);
/*@}*/

extern p3m_struct p3m;
/*----------------------------------------------------------*/

#endif
