#ifndef ELC_H
#define ELC_H

typedef struct {
  double maxPWerror;
  double far_cut, far_cut2;
  double minimal_dist;
  int far_calculated;
} ELC_struct;
extern ELC_struct elc_params;

/** set parameters for ELC */
int set_elc_params(Tcl_Interp *interp, double maxPWerror, double min_dist, double far_cut);

void ELC_add_force();

double ELC_energy();

void ELC_init();

void ELC_on_resort_particles();

#endif
