#ifndef MYP3M_H
#define MYP3M_H


/*----------------------------------------------------------*/

typedef struct {
  double alpha;     /** Ewald splitting parameter */
  int    mesh[3];      /** number of mesh points per coordinate direction */
  int    P;         /** charge assignment order */
  double rcut;      /** cutoff radius for real space electrostatics */
  double prefactor; /** Coulomb prefactor, e.g. Bjerrum*Temp */
  double epsilon;   /** epsilon of the "surrounding dielectric" */
} p3m_struct;

extern double skin; 
extern void   setp3m();
extern void   P3M_init();
extern void   compute_measures();
extern void   P3M_perform();
extern void   P3M_exit();
extern void interpolate_charge_assignment_function(void);

#endif
