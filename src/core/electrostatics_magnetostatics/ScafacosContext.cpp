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
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/range/numeric.hpp>

#include <cassert>

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
    auto const G = Utils::Vector<Utils::Vector3d, 3>{
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

} // namespace Scafacos
#endif /* SCAFACOS */
