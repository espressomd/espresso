#ifndef OBSERVABLES_PROFIELS_HPP
#define OBSERVABLES_PROFIELS_HPP

#include "utils/List.hpp"

typedef struct {
  IntList *id_list;
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_z;
  double max_z;
  int n_x_bins;
  int n_y_bins;
  int n_z_bins;
  void *container;
} profile_data;

typedef struct {
  IntList *id_list;
  double min_r;
  double max_r;
  double min_phi;
  double max_phi;
  double min_z;
  double max_z;
  double center[3];
  double axis[3];
  int n_phi_bins;
  int n_r_bins;
  int n_z_bins;
  void *container;
} radial_profile_data;

#endif
