#ifndef MMM1D_H
#define MMM1D_H
typedef struct {
  double far_switch_radius_2;
  int    bessel_cutoff;
  double maxPWerror;
  double prefactor;
  double bjerrum;
} MMM1D_struct;
extern MMM1D_struct mmm1d;

void set_mmm1d_params(double bjerrum, double switch_rad,
		      int bessel_cutoff, double maxPWerror);

void MMM1D_recalcTables();

void MMM1D_calc_forces();

void MMM1D_init();
#endif
