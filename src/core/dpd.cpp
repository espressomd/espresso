/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
     Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file dpd.cpp
    Implementation of \ref dpd.hpp "dpd.hpp"
 */
#include "dpd.hpp"

#ifdef DPD

#include "integrate.hpp"
#include "random.hpp"
#include "thermostat.hpp"

#ifdef LEES_EDWARDS
/** Chatterjee 2007 proposes that for DPD with Lees Edwards BCs,
 *  it is better not to count interactions with Ghost particles. */
static bool le_chatterjee_test_pair(Particle *p1, Particle *p2) {
  /* cannot assume that we have a flag set to mark ghost particles,
   * but can assume (for LE) that non-ghost
   * particle y-coords are always imaged inside the box */
  if (p2->l.ghost or p1->l.ghost)
    return false;
  else
    return true;
}
#endif

void dpd_heat_up() {
  double pref_scale = sqrt(3);
  dpd_update_params(pref_scale);
}

void dpd_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  dpd_update_params(pref_scale);
}

int dpd_set_params(int part_type_a, int part_type_b, double gamma, double r_c,
                   int wf, double tgamma, double tr_c, int twf) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  data->dpd_gamma = gamma;
  data->dpd_r_cut = r_c;
  data->dpd_wf = wf;
  data->dpd_pref2 = sqrt(24.0 * temperature * gamma / time_step);
  data->dpd_tgamma = tgamma;
  data->dpd_tr_cut = tr_c;
  data->dpd_twf = twf;
  data->dpd_pref4 = sqrt(24.0 * temperature * tgamma / time_step);

  /* Only make active if the DPD thermostat is
     activated, otherwise it will by activated
     by dpd_init() on thermostat change.
  */
  if (thermo_switch & THERMO_DPD) {
    data->dpd_pref1 = gamma / time_step;
    data->dpd_pref3 = tgamma / time_step;
  } else {
    data->dpd_pref1 = 0.0;
    data->dpd_pref3 = 0.0;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

void dpd_init() {
  for (int type_a = 0; type_a < n_particle_types; type_a++) {
    for (int type_b = 0; type_b < n_particle_types; type_b++) {
      auto data = get_ia_param(type_a, type_b);
      if ((data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0)) {
        data->dpd_pref1 = data->dpd_gamma / time_step;
        data->dpd_pref2 =
            sqrt(24.0 * temperature * data->dpd_gamma / time_step);
        data->dpd_pref3 = data->dpd_tgamma / time_step;
        data->dpd_pref4 =
            sqrt(24.0 * temperature * data->dpd_tgamma / time_step);
      }
    }
  }
}

void dpd_switch_off(void) {
  for (int type_a = 0; type_a < n_particle_types; type_a++) {
    for (int type_b = 0; type_b < n_particle_types; type_b++) {
      auto data = get_ia_param(type_a, type_b);
      data->dpd_pref1 = data->dpd_pref3 = 0.0;
    }
  }
}

void dpd_update_params(double pref_scale) {
  int type_a, type_b;
  IA_parameters *data;

  for (type_a = 0; type_a < n_particle_types; type_a++) {
    for (type_b = 0; type_b < n_particle_types; type_b++) {
      data = get_ia_param(type_a, type_b);
      if ((data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0)) {
        data->dpd_pref2 *= pref_scale;
        data->dpd_pref4 *= pref_scale;
      }
    }
  }
}

void add_dpd_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
                        double d[3], double dist, double dist2) {
  int j;
  // velocity difference between p1 and p2
  double vel12_dot_d12 = 0.0;
  // weighting functions for friction and random force
  double omega, omega2; // omega = w_R/dist
  double friction, noise;
  // Projection martix
  int i;
  double P_times_dist_sqr[3][3] = {{dist2, 0, 0}, {0, dist2, 0}, {0, 0, dist2}},
         noise_vec[3];
  double f_D[3], f_R[3];
  double tmp;

#ifdef LEES_EDWARDS
  if (le_chatterjee_test_pair(p1, p2))
    return;
#endif

  auto const dist_inv = 1.0 / dist;

  if ((dist < ia_params->dpd_r_cut) && (ia_params->dpd_pref1 > 0.0)) {
    if (ia_params->dpd_wf == 0) {
      omega = dist_inv;
    } else {
      omega = dist_inv - 1.0 / ia_params->dpd_r_cut;
    }

    omega2 = SQR(omega);
    // DPD part
    // friction force prefactor
    for (j = 0; j < 3; j++)
      vel12_dot_d12 += (p1->m.v[j] - p2->m.v[j]) * d[j];
    friction = ia_params->dpd_pref1 * omega2 * vel12_dot_d12;
    // random force prefactor
    noise = ia_params->dpd_pref2 * omega * (d_random() - 0.5);
    for (j = 0; j < 3; j++) {
      p1->f.f[j] += (tmp = (noise - friction) * d[j]);
      p2->f.f[j] -= tmp;
    }
  }
  // DPD2 part
  if ((dist < ia_params->dpd_tr_cut) && (ia_params->dpd_pref3 > 0.0)) {
    if (ia_params->dpd_twf == 0) {
      omega = dist_inv;
    } else {
      omega = dist_inv - 1.0 / ia_params->dpd_tr_cut;
    }

    omega2 = SQR(omega);

    for (i = 0; i < 3; i++) {
      // noise vector
      noise_vec[i] = d_random() - 0.5;
      // Projection Matrix
      for (j = 0; j < 3; j++) {
        P_times_dist_sqr[i][j] -= d[i] * d[j];
      }
    }
    for (i = 0; i < 3; i++) {
      // Damping force
      f_D[i] = 0;
      // Random force
      f_R[i] = 0;
      for (j = 0; j < 3; j++) {
        f_D[i] += P_times_dist_sqr[i][j] * (p1->m.v[j] - p2->m.v[j]);
        f_R[i] += P_times_dist_sqr[i][j] * noise_vec[j];
      }
      // NOTE: velocity are scaled with time_step
      f_D[i] *= ia_params->dpd_pref3 * omega2;
      // NOTE: noise force scales with 1/sqrt(time_step
      f_R[i] *= ia_params->dpd_pref4 * omega * dist_inv;
    }
    for (j = 0; j < 3; j++) {
      tmp = f_R[j] - f_D[j];
      p1->f.f[j] += tmp;
      p2->f.f[j] -= tmp;
    }
  }
}

#endif
