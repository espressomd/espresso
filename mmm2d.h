#ifndef MMM2D_H
#define MMM2D_H

typedef struct {
  double maxPWerror;
  double far_cut, far_cut2;
  int far_calculated;
} MMM2D_struct;
extern MMM2D_struct mmm2d_params;

/** set parameters for MMM2D */
int set_mmm2d_params(Tcl_Interp *interp, double maxPWerror, double far_cut);

/** the actual long range force calculation */
void MMM2D_add_far_force();

/** the actual long range energy calculation */
double MMM2D_far_energy();

/** pairwise calculated parts of MMM2D force (near neighbors) */
void add_mmm2d_coulomb_pair_force(Particle *p1, Particle *p2,
				  double dv[3], double d2, double d);

/** pairwise calculated parts of MMM2D force (near neighbors) */
double mmm2d_coulomb_pair_energy(Particle *p1, Particle *p2,
				 double dv[3], double d2, double d);

void MMM2D_init();

void MMM2D_on_resort_particles();

#endif
