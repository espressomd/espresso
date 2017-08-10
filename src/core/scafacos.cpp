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

#include "scafacos.hpp"

#if defined(SCAFACOS)

#include <vector>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "scafacos/Scafacos.hpp"
#include "tuning.hpp"

/** This file contains the c-like interface for Scafacos */

namespace Scafacos {

/** Get available scafacos methods */
std::list<std::string> available_methods() {
  return Scafacos::available_methods();
}

/** Encapsulation for the particle data needed by scafacos */
struct ScafacosData {
  int update_particle_data();
  void update_particle_forces() const;
  /** Inputs */
  std::vector<double> charges, positions, dipoles;
  /** Outputs */
  std::vector<double> fields, potentials;
};

static Scafacos *scafacos = 0;
static ScafacosData particles;

/** \brief Collect particle data in continous arrays as required by fcs */
int ScafacosData::update_particle_data() {
  positions.clear();

  if (!dipolar()) {
    charges.clear();
  } else {
#ifdef SCAFACOS_DIPOLES
    dipoles.clear();
#endif
  }

  for (auto const &p : local_cells.particles()) {
    positions.push_back(p.r.p[0]);
    positions.push_back(p.r.p[1]);
    positions.push_back(p.r.p[2]);
    if (!dipolar()) {
      charges.push_back(p.p.q);
    } else {
#ifdef SCAFACOS_DIPOLES
      dipoles.push_back(p.r.dip[0]);
      dipoles.push_back(p.r.dip[1]);
      dipoles.push_back(p.r.dip[2]);
#endif
    }
  }

  return positions.size() / 3;
}

/** \brief Write forces back to particles */

void ScafacosData::update_particle_forces() const {
  int it = 0;

  for (auto &p : local_cells.particles()) {
    if (!dipolar()) {
      p.f.f[0] += coulomb.prefactor * p.p.q * fields[it++];
      p.f.f[1] += coulomb.prefactor * p.p.q * fields[it++];
      p.f.f[2] += coulomb.prefactor * p.p.q * fields[it++];
    } else {
#ifdef SCAFACOS_DIPOLES
      // Indices
      // 3 "potential" values per particles (see below)
      int it_p = 3 * it;
      // 6 "field" values per particles (see below)
      int it_f = 6 * it;

      // The scafacos term "potential" here in fact refers to the magnetic
      // field
      // So, the torques are given by m \times B
      double t[3];
      utils::cross_product(p.r.dip, &(potentials[it_p]), t);
      // The force is given by G m, where G is a matrix
      // which comes from teh "fields" output of scafacos like this
      // 0 1 2
      // 1 3 4
      // 2 4 5
      // where the numbers refer to indices in the "field" output from
      // scafacos
      double f[3];
      f[0] = fields[it_f + 0] * p.r.dip[0] + fields[it_f + 1] * p.r.dip[1] +
             fields[it_f + 2] * p.r.dip[2];
      f[1] = fields[it_f + 1] * p.r.dip[0] + fields[it_f + 3] * p.r.dip[1] +
             fields[it_f + 4] * p.r.dip[2];
      f[2] = fields[it_f + 2] * p.r.dip[0] + fields[it_f + 4] * p.r.dip[1] +
             fields[it_f + 5] * p.r.dip[2];

      // Add to particles
      for (int j = 0; j < 3; j++) {
        p.f.f[j] += coulomb.Dprefactor * f[j];
        p.f.torque[j] += coulomb.Dprefactor * t[j];
      }
      it++;
#endif
    }
  }

  /** Check that the particle number did not change */
  if (!dipolar()) {
    assert(it == fields.size());
  } else {
    assert(it = positions.size() / 3);
  }
}

void add_pair_force(Particle *p1, Particle *p2, double *d, double dist,
                    double *force) {
  if (dist > get_r_cut())
    return;

  assert(scafacos);
  const double field = scafacos->pair_force(dist);
  const double fak = p2->p.q * p1->p.q * field * coulomb.prefactor / dist;

  for (int i = 0; i < 3; i++) {
    p1->f.f[i] -= fak * d[i];
    p2->f.f[i] += fak * d[i];
  }
}

double pair_energy(Particle *p1, Particle *p2, double dist) {
  if (dist <= get_r_cut())
    return coulomb.prefactor * p1->p.q * p2->p.q * scafacos->pair_energy(dist);
  else
    return 0.;
}

void add_long_range_force() {
  particles.update_particle_data();

  if (scafacos) {
    if (!dipolar()) {
      scafacos->run(particles.charges, particles.positions, particles.fields,
                    particles.potentials);
    } else {
#ifdef SCAFACOS_DIPOLES
      scafacos->run_dipolar(particles.dipoles, particles.positions,
                            particles.fields, particles.potentials);
#endif
    }
  } else
    runtimeError("Scafacos internal error.");

  particles.update_particle_forces();
}

double long_range_energy() {
  if (scafacos) {
    particles.update_particle_data();
    if (!dipolar()) {
      scafacos->run(particles.charges, particles.positions, particles.fields,
                    particles.potentials);
      return 0.5 * coulomb.prefactor *
             std::inner_product(particles.charges.begin(),
                                particles.charges.end(),
                                particles.potentials.begin(), 0.0);
    } else {
#ifdef SCAFACOS_DIPOLES
      scafacos->run_dipolar(particles.dipoles, particles.positions, particles.fields,particles.potentials);
      return -0.5* coulomb.Dprefactor * std::inner_product(particles.dipoles.begin(),particles.dipoles.end(),particles.potentials.begin(),0.0);
#endif
    }
  }

  return 0.0;
}

/** Determine runtime for a specific cutoff */
double time_r_cut(double r_cut) {
  double t;

  /** Set cutoff to time */
  scafacos->set_r_cut(r_cut);

  /** Tune other parameters */
  scafacos->tune(particles.charges, particles.positions);

  on_coulomb_change();

  return time_force_calc(10);
}

/** Determine the optimal cutoff by bisection */
void tune_r_cut() {
  const double tune_limit = 1e-3;
  /** scafacos p3m and ewald do not accept r_cut 0 for no good reason */
  double r_min = 1.0;
  double r_max = std::min(min_local_box_l, min_box_l / 2.0) - skin;
  double t_min = 0;
  double t_max = std::numeric_limits<double>::max();

  /** Run bisection */
  while (std::fabs(r_min - r_max) > tune_limit) {
    const double dr = 0.5 * (r_max - r_min);
    const double t_mid = time_r_cut(r_min + dr);
    t_min = time_r_cut(r_min);
    t_max = time_r_cut(r_max);

    if (t_min <= 0.0) {
      r_min += tune_limit;
      break;
    }

    if (t_max <= 0.0) {
      r_min -= tune_limit;
      break;
    }

    if (t_mid > t_min) {
      r_max = r_min += dr;
    } else {
      r_min += dr;
    }
  }
}

void tune() {
  particles.update_particle_data();

  /** Check whether we have to do a bisection for the short range cutoff */
  /** Check if there is a user supplied cutoff */
  if ((scafacos->has_near) && (scafacos->r_cut() <= 0.0)) {
    tune_r_cut();
  } else {
    scafacos->tune(particles.charges, particles.positions);
  }
}

static void set_params_safe(const std::string &method, const std::string &params, bool dipolar_ia) {
  if(scafacos) {
    delete scafacos;
    scafacos = 0;
  }

  scafacos = new Scafacos(method, comm_cart, params);

  int per[3] = { PERIODIC(0) != 0, PERIODIC(1) != 0, PERIODIC(2) != 0 };

  scafacos->set_dipolar(dipolar_ia);
  scafacos->set_common_parameters(box_l, per, n_part);

  on_coulomb_change();

  if (!dipolar_ia) {
    tune();

    on_coulomb_change();
  }
}

/** Bend result from scafacos back to original format */
std::string get_method_and_parameters() {
  if(!scafacos) {
    return std::string();
  }

  std::string p = scafacos->get_method()+" "+scafacos->get_parameters();

  std::replace(p.begin(), p.end(), ',', ' ');

  return p;
}

double get_r_cut() {
  if(scafacos) {
    if (!scafacos->has_near) return 0;
    return scafacos->r_cut();
  }
  return 0.0;
}

void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar_ia) {
  mpi_call(mpi_scafacos_set_parameters_slave, method.size(), params.size());

  /** This requires C++11, otherwise this is undefined because std::string was
   * not required to have conitnuous memory before. */
  /* const_cast is ok, this code runs only on rank 0 where the mpi call does not
   * modify the buffer */
  MPI_Bcast(const_cast<char *>(&(*method.begin())), method.size(), MPI_CHAR, 0,
            comm_cart);
  MPI_Bcast(const_cast<char *>(&(*params.begin())), params.size(), MPI_CHAR, 0,
            comm_cart);

#ifdef SCAFACOS_DIPOLES
  bool d = dipolar_ia;
  MPI_Bcast(&d, sizeof(bool), MPI_CHAR, 0, comm_cart);
#endif

  set_params_safe(method, params, dipolar_ia);
#ifdef SCAFACOS_DIPOLES
  set_dipolar(d);
#endif
}

bool dipolar() {
  if (scafacos)
    return scafacos->dipolar();
  runtimeErrorMsg() << "Scafacos not initialized";
}

void set_dipolar(bool d) {
  if (scafacos) {
    scafacos->set_dipolar(d);
  } else
    runtimeErrorMsg() << "Scafacos not initialized.";
}

void free_handle() {

  if (this_node==0) 
    mpi_call(mpi_scafacos_free_slave, 0,0);
  if(scafacos) {
delete scafacos;
    scafacos = 0;
  }
}

void on_boxl_change() {
// If scafacos is not active, do nothing
if (!scafacos) return;

// Get current parameters and re_initialize
std::string params=scafacos->get_parameters();
std::string method=scafacos->get_method();
bool dip=scafacos->dipolar();

// Delete existing scafacos instance
free_handle();

// And make a new one
set_parameters(method,params,dip);
}

} // namespace scafacos
#endif /* SCAFACOS */

void mpi_scafacos_set_parameters_slave(int n_method, int n_params) {
#if defined(SCAFACOS)
  using namespace Scafacos;
  std::string method;
  std::string params;

  method.resize(n_method);
  params.resize(n_params);

  /** This requires C++11, otherwise this is undefined because std::string was
   * not required to have conitnuous memory before. */
  MPI_Bcast(&(*method.begin()), n_method, MPI_CHAR, 0, comm_cart);
  MPI_Bcast(&(*params.begin()), n_params, MPI_CHAR, 0, comm_cart);
  bool dip = false;
#ifdef SCAFACOS_DIPOLES
  MPI_Bcast(&dip, sizeof(bool), MPI_CHAR, 0, comm_cart);
#endif

  set_params_safe(method, params, dip);
#ifdef SCAFACOS_DIPOLES
  set_dipolar(dip);
#endif
#endif /* SCAFACOS */
}

void mpi_scafacos_free_slave(int a, int b) {
  #if defined(SCAFACOS) 
  using namespace Scafacos;
  free_handle();
  #endif
}
