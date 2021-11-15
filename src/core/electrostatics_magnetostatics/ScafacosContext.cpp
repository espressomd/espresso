/*
 * Copyright (C) 2010-2020 The ESPResSo project
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

#include "config.hpp"

#if defined(SCAFACOS)

#include "electrostatics_magnetostatics/ScafacosContext.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "electrostatics_magnetostatics/scafacos.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/range/algorithm/min_element.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <limits>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

namespace Scafacos {

std::string ScafacosContext::get_method_and_parameters() {
  std::string representation = get_method() + " " + get_parameters();
  std::replace(representation.begin(), representation.end(), ',', ' ');
  return representation;
}

void ScafacosContext::update_system_params() {
  int per[3] = {box_geo.periodic(0) != 0, box_geo.periodic(1) != 0,
                box_geo.periodic(2) != 0};

  auto const n_part = boost::mpi::all_reduce(
      comm_cart, cell_structure.local_particles().size(), std::plus<>());
  set_common_parameters(box_geo.length().data(), per, n_part);
}

void ScafacosContextCoulomb::update_particle_data() {
  positions.clear();
  charges.clear();

  for (auto const &p : cell_structure.local_particles()) {
    auto const pos = folded_position(p.r.p, box_geo);
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    charges.push_back(p.p.q);
  }
}

#ifdef SCAFACOS_DIPOLES
void ScafacosContextDipoles::update_particle_data() {
  positions.clear();
  dipoles.clear();

  for (auto const &p : cell_structure.local_particles()) {
    auto const pos = folded_position(p.r.p, box_geo);
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    auto const dip = p.calc_dip();
    dipoles.push_back(dip[0]);
    dipoles.push_back(dip[1]);
    dipoles.push_back(dip[2]);
  }
}
#endif

void ScafacosContextCoulomb::update_particle_forces() const {
  if (positions.empty())
    return;

  int it = 0;
  for (auto &p : cell_structure.local_particles()) {
    p.f.f += coulomb.prefactor * p.p.q *
             Utils::Vector3d(Utils::Span<const double>(&(fields[it]), 3));
    it += 3;
  }

  /* Check that the particle number did not change */
  assert(it == fields.size());
}

#ifdef SCAFACOS_DIPOLES
void ScafacosContextDipoles::update_particle_forces() const {
  if (positions.empty())
    return;

  int it = 0;
  for (auto &p : cell_structure.local_particles()) {
    // Indices
    // 3 "potential" values per particles (see below)
    int const it_p = 3 * it;
    // 6 "field" values per particles (see below)
    int const it_f = 6 * it;

    // The scafacos term "potential" here in fact refers to the magnetic
    // field. So, the torques are given by m \times B
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
    auto const G = Utils::Matrix<double, 3, 3>{
        {fields[it_f + 0], fields[it_f + 1], fields[it_f + 2]},
        {fields[it_f + 1], fields[it_f + 3], fields[it_f + 4]},
        {fields[it_f + 2], fields[it_f + 4], fields[it_f + 5]}};
    auto const f = G * dip;

    // Add to particles
    p.f.f += dipole.prefactor * f;
    p.f.torque += dipole.prefactor * t;
    it++;
  }

  /* Check that the particle number did not change */
  assert(it == positions.size() / 3);
}
#endif

double ScafacosContextCoulomb::long_range_energy() {
  update_particle_data();
  run(charges, positions, fields, potentials);
  return 0.5 * coulomb.prefactor *
         boost::inner_product(charges, potentials, 0.0);
}

#ifdef SCAFACOS_DIPOLES
double ScafacosContextDipoles::long_range_energy() {
  update_particle_data();
  run_dipolar(dipoles, positions, fields, potentials);
  return -0.5 * dipole.prefactor *
         boost::inner_product(dipoles, potentials, 0.0);
}
#endif

static void set_r_cut_and_tune_local(double r_cut) {
  set_r_cut_and_tune(r_cut);
}

REGISTER_CALLBACK(set_r_cut_and_tune_local)

/** Determine runtime for a specific cutoff */
static double time_r_cut(double r_cut) {
  /* Set cutoff to time */
  mpi_call_all(set_r_cut_and_tune_local, r_cut);

  return time_force_calc(10);
}

/** Determine the optimal cutoff by bisection */
static void tune_r_cut() {
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

void ScafacosContextCoulomb::tune() {
  update_particle_data();

  // Check whether we have to do a bisection for the short range cutoff
  // Check if there is a user-supplied cutoff
  if (has_near && r_cut() <= 0.0) {
    // run tuning on the head node
    if (this_node == 0) {
      tune_r_cut();
    }
  } else {
    // ESPResSo is not affected by a short range cutoff -> tune in parallel
    Scafacos::tune(charges, positions);
  }
}

void ScafacosContextCoulomb::set_r_cut_and_tune(double r_cut) {
  update_particle_data();
  set_r_cut(r_cut);
  Scafacos::tune(charges, positions);
}

} // namespace Scafacos
#endif /* SCAFACOS */
