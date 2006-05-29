/*****************************************************************************
 *                   basic types for machine independence                  
 ****************************************************************************/

#ifndef LB_H
#define LB_H

#include <tcl.h>
#include "utils.h"
#include "grid.h"
#include "global.h"

#ifdef LB

typedef int                LBT_BOOL; 

///
struct LB_structure {
  /** Number of gridpoints in each dimension  */
  int     gridpoints;   
  /** number density [1/a^3] */
  double  rho;          /* number density [1/a^3] */
  /** kinematic viscosity [a^2/tau] */
  double  viscos;      /* kinematic viscosity [a^2/tau] */
  /** lb frinction coefficient. Usually chosen to be 20 */
  double  lbfriction;
  /** grid distance [LJ units] */
  double  agrid;        /* grid distance [LJ units] */
  /** time step to propagate from one gridpoint to next [LJ] also used for Verlet step (preleminary */
  double  tau;         /* time step to propagate from one gridpoint to next [LJ] also used for Verlet step (preleminary */

  double  lblambda;      /* relaxation time [tau] for non-eq distribution, not independent from viscos */
  double  c_sound_sq;  /* (sound velocity)^2 [a^2/tau^2]  */
  double  current_type;     /* current density [1/a^2/tau] */
  char    currentpar[200];
  LBT_BOOL  boundary;    /* boundaries yes/no */
  int     chainlen;        
          
};
typedef struct LB_structure LB_structure;

extern LB_structure lbpar; 

/** If non-zero, the lb populations will not be updated during force recalc. */
extern int transfer_momentum;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

void calc_fluid_chain_interaction(int iamghost, int nMono);
/**
  calculation of hydrodynamic interaction. fluid and monomer velocities
  are coupled by a friction coefficient
  
  input: friction, position of monomers (xp), hydrodynamic fields (rho,j,n)
  stochastic force given to monomers (forces)
	
  output: hydrodynamic force is added to variable 'force' 
  and via momentum conservation: shifting of population n
  this works only in three dimensions
*/

int  lb_callback(Tcl_Interp *interp, void *_data);

void calc_lbforce();
/**
   This function is used to calculate fluid particle interaction
*/
 
void LB_propagate();
/**
   Main function for propagating the LB part of the code. It consists of two 
   parts: collion and propagation.
*/
int  lb_parser(Tcl_Interp * interp, int argc, char ** argv);
int  printlbIAToResult(Tcl_Interp *interp);
void thermo_init_lb();

/*@}*/

#endif /* LB */

#endif /* LB_H */







