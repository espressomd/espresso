/*------------------------------------------------------------

HEADER-FILE fuer p3m_v2.c

Letzte Aenderung: 18.11.98

Markus Deserno

------------------------------------------------------------*/

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

/*----------------------------------------------------------*/

typedef struct {
  double alpha;     /* Ewald splitting parameter */
  int    mesh;      /* number of mesh points per coordinate direction */
  int    P;         /* charge assignment order */
  double rcut;      /* cutoff radius for real space electrostatics */
  double prefactor; /* Coulomb prefactor, e.g. Bjerrum*Temp */
  double epsilon;   /* epsilon of the "surrounding dielectric" */
}p3m_struct;

/*----------------------------------------------------------*/

extern void   P3M_init(double *, int, p3m_struct *p3m);

extern void   P3M_exit(void);

extern double P3M_dipole(int force);

extern void   P3M_perform(int force, double *E_Coulomb_P3M);

extern double P3M_perform_subset(int, int);

extern p3m_struct p3m;
/*----------------------------------------------------------*/

#endif
