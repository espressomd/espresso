/*****************************************************************************/
/*                     basic types for machine independence                  */
/*****************************************************************************/

#ifdef LB

#ifndef LB_H
#define LB_H
typedef int                LBT_BOOL; 

typedef struct {
  int     gridpoints;   /* # of gridpoints in each dimension  */
  double  rho;          /* number density [1/a^3] */
  double  viscos;      /* kinematic viscosity [a^2/tau] */
  double  lbfriction;
  double  agrid;        /* grid distance [LJ units] */
  double  tau;         /* time step to propagate from one gridpoint to next [LJ] also used for Verlet step (preleminary */

  double  lblambda;      /* relaxation time [tau] for non-eq distribution, not independent from viscos */
  double  c_sound_sq;  /* (sound velocity)^2 [a^2/tau^2]  */
  double  current_type;     /* current density [1/a^2/tau] */
  char    currentpar[200];
  int     ranseed;     /* seed for random number generator   */
  LBT_BOOL  boundary;    /* boundaries yes/no */
  int     chainlen;        
          
} LB_structure;

extern LB_structure lbpar; 

/** If non-zero, the lb populations will no be updated during force recalc. */
extern int transfer_momentum;

void LB_Run();
void calc_fluid_chain_interaction(int iamghost);
int  lb_callback(Tcl_Interp *interp, void *_data);

void calc_lbforce();
void LB_propagate();
int  lb_parser(Tcl_Interp * interp, int argc, char ** argv);
int  printlbIAToResult(Tcl_Interp *interp);
void thermo_init_lb();

#endif /* LB_H */

#endif /* LB */







