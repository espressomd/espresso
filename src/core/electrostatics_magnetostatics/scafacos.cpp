/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** @file
 *  Provide a C-like interface for ScaFaCoS.
 */

#include "electrostatics_magnetostatics/scafacos.hpp"

#if defined(SCAFACOS)

#include <boost/mpi/collectives.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "Scafacos.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "tuning.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#if defined(SCAFACOS_DIPOLES) && !defined(FCS_ENABLE_DIPOLES)
#error                                                                         \
    "SCAFACOS_DIPOLES requires dipoles support in scafacos library (FCS_ENABLE_DIPOLES)."
#endif

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

static Scafacos *scafacos = nullptr;
static ScafacosData particles;

/** \brief Collect particle data in continuous arrays as required by fcs */
int ScafacosData::update_particle_data() {
  positions.clear();

  if (!dipolar()) {
    charges.clear();
  } else {
#ifdef SCAFACOS_DIPOLES
    dipoles.clear();
#endif
  }

  for (auto const &p : cell_structure.local_particles()) {
    auto const pos = folded_position(p.r.p, box_geo);
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    if (!dipolar()) {
      charges.push_back(p.p.q);
    } else {
#ifdef SCAFACOS_DIPOLES
      auto const dip = p.calc_dip();
      dipoles.push_back(dip[0]);
      dipoles.push_back(dip[1]);
      dipoles.push_back(dip[2]);
#endif
    }
  }

  return static_cast<int>(positions.size() / 3);
}

/** \brief Write forces back to particles */
void ScafacosData::update_particle_forces() const {
  int it = 0;
  if (positions.empty())
    return;

  for (auto &p : cell_structure.local_particles()) {
    if (!dipolar()) {
      p.f.f += coulomb.prefactor * p.p.q *
               Utils::Vector3d(Utils::Span<const double>(&(fields[it]), 3));
      it += 3;
    } else {
#ifdef SCAFACOS_DIPOLES
      // Indices
      // 3 "potential" values per particles (see below)
      int const it_p = 3 * it;
      // 6 "field" values per particles (see below)
      int const it_f = 6 * it;

      // The scafacos term "potential" here in fact refers to the magnetic
      // field
      // So, the torques are given by m \times B
      auto const dip = p.calc_dip();
      auto const t = vector_product(
          dip,
          Utils::Vector3d(Utils::Span<const double>(&(potentials[it_p]), 3)));
      // The force is given by G m, where G is a matrix
      // which comes from the "fields" output of scafacos like this
      // 0 1 2
      // 1 3 4
      // 2 4 5
      // where the numbers refer to indices in the "field" output from scafacos
      auto const G = Utils::Vector<Utils::Vector3d, 3>{
          {fields[it_f + 0], fields[it_f + 1], fields[it_f + 2]},
          {fields[it_f + 1], fields[it_f + 3], fields[it_f + 4]},
          {fields[it_f + 2], fields[it_f + 4], fields[it_f + 5]}};
      auto const f = G * dip;

      // Add to particles
      p.f.f += dipole.prefactor * f;
      p.f.torque += dipole.prefactor * t;
      it++;
#endif
    }
  }

  /* Check that the particle number did not change */
  if (!dipolar()) {
    assert(it == fields.size());
  } else {
    assert(it == positions.size() / 3);
  }
}

void add_pair_force(double q1q2, const double *d, double dist, double *force) {
  if (dist > get_r_cut())
    return;

  assert(scafacos);
  const double field = scafacos->pair_force(dist);
  const double fak = q1q2 * field / dist;

  for (int i = 0; i < 3; i++) {
    force[i] -= fak * d[i];
  }
}

double pair_energy(double q1q2, double dist) {
  if (dist <= get_r_cut())
    return q1q2 * scafacos->pair_energy(dist);
  return 0.;
}

// Issues a runtime error if positions are outside the box domain
// This is needed, because the scafacos grid sort produces an mpi deadlock
// otherwise
// Returns true if calculations can continue.
bool check_position_validity(const std::vector<double> &pos) {
  assert(pos.size() % 3 == 0);
  for (int i = 0; i < pos.size() / 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (pos[3 * i + j] < 0 || pos[3 * i + j] > box_geo.length()[j]) {
        // Throwing exception rather than runtime error, because continuing will
        // result in mpi deadlock
        throw std::runtime_error("Particle position outside the box domain not "
                                 "allowed for scafacos-based methods.");
      }
    }
  }
  return true;
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
    throw std::runtime_error(
        "Scafacos internal error. Instance pointer is not valid.");

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
    }
#ifdef SCAFACOS_DIPOLES
    scafacos->run_dipolar(particles.dipoles, particles.positions,
                          particles.fields, particles.potentials);
    return -0.5 * dipole.prefactor *
           std::inner_product(particles.dipoles.begin(),
                              particles.dipoles.end(),
                              particles.potentials.begin(), 0.0);
#endif
  }

  return 0.0;
}
void set_r_cut_and_tune_local(double r_cut) {
  particles.update_particle_data();

  scafacos->set_r_cut(r_cut);
  scafacos->tune(particles.charges, particles.positions);
}

REGISTER_CALLBACK(set_r_cut_and_tune_local)

/** Determine runtime for a specific cutoff */
double time_r_cut(double r_cut) {
  /* Set cutoff to time */
  mpi_call(set_r_cut_and_tune_local, r_cut);
  set_r_cut_and_tune_local(r_cut);

  return time_force_calc(10);
}

/** Determine the optimal cutoff by bisection */
void tune_r_cut() {
  const double tune_limit = 1e-3;

  auto const min_box_l = *boost::min_element(box_geo.length());
  auto const min_local_box_l = *boost::min_element(local_geo.length());

  /* scafacos p3m and Ewald do not accept r_cut 0 for no good reason */
  double r_min = 1.0;
  double r_max = std::min(min_local_box_l, min_box_l / 2.0) - skin;
  double t_min = 0;
  double t_max = std::numeric_limits<double>::max();

  /* Run bisection */
  while (std::fabs(r_min - r_max) > tune_limit) {
    const double dr = 0.5 * (r_max - r_min);
    const double t_mid = time_r_cut(r_min + dr);
    t_min = time_r_cut(r_min);
    t_max = time_r_cut(r_max);

    if (t_min <= 0.0 or t_max <= 0.0) {
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

  /* Check whether we have to do a bisection for the short range cutoff */
  /* Check if there is a user supplied cutoff */
  if ((scafacos->has_near) && (scafacos->r_cut() <= 0.0)) {
    // Tuning of r_cut needs to run on the master node because it relies on
    // master-slave mode communication
    if (this_node == 0) {
      tune_r_cut();
    } else {
      return; // Tune on the master node will issue mpi calls
    }
  } else {
    // ESPResSo is not affected by a short range cutoff. Tune in parallel
    scafacos->tune(particles.charges, particles.positions);
  }
}

static void set_params_safe(const std::string &method,
                            const std::string &params, bool dipolar_ia,
                            int n_part) {
  if (scafacos) {
    delete scafacos;
    scafacos = nullptr;
  }

  scafacos = new Scafacos(method, comm_cart, params);

  int per[3] = {box_geo.periodic(0) != 0, box_geo.periodic(1) != 0,
                box_geo.periodic(2) != 0};

  scafacos->set_dipolar(dipolar_ia);
#ifdef DIPOLES
  if (dipolar_ia) {
    dipole.method = DIPOLAR_SCAFACOS;
  }
#endif
#ifdef ELECTROSTATICS
  if (!dipolar_ia) {
    coulomb.method = COULOMB_SCAFACOS;
  }
#endif
  scafacos->set_common_parameters(box_geo.length().data(), per, n_part);

  on_coulomb_change();

  if (!dipolar_ia) {
    tune();
  }
}

REGISTER_CALLBACK(set_params_safe)

/** Bend result from scafacos back to original format */
std::string get_method_and_parameters() {
  if (!scafacos) {
    return std::string();
  }

  std::string p = scafacos->get_method() + " " + scafacos->get_parameters();

  std::replace(p.begin(), p.end(), ',', ' ');

  return p;
}

double get_r_cut() {
  if (scafacos) {
    if (!scafacos->has_near)
      return 0;
    return scafacos->r_cut();
  }
  return 0.0;
}

void set_parameters(const std::string &method, const std::string &params,
                    bool dipolar_ia) {
  mpi_call_all(set_params_safe, method, params, dipolar_ia, get_n_part());
}

bool dipolar() {
  if (scafacos)
    return scafacos->dipolar();
  throw std::runtime_error("Scafacos not initialized");
}

void set_dipolar(bool d) {
  if (scafacos) {
    scafacos->set_dipolar(d);
  } else
    runtimeErrorMsg() << "Scafacos not initialized.";
}

void free_handle() {
  if (this_node == 0)
    mpi_call(free_handle);
  if (scafacos) {
    delete scafacos;
    scafacos = nullptr;
  }
}

REGISTER_CALLBACK(free_handle)

void update_system_params() {
  // If scafacos is not active, do nothing
  if (!scafacos) {
    throw std::runtime_error("Scafacos object not there");
  }

  int per[3] = {box_geo.periodic(0) != 0, box_geo.periodic(1) != 0,
                box_geo.periodic(2) != 0};

  auto const n_part = boost::mpi::all_reduce(comm_cart, cells_get_n_particles(),
                                             std::plus<int>());
  scafacos->set_common_parameters(box_geo.length().data(), per, n_part);
}

} // namespace Scafacos
#endif /* SCAFACOS */
