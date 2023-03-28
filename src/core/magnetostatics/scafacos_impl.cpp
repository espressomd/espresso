/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef SCAFACOS_DIPOLES

#include "magnetostatics/scafacos.hpp"
#include "magnetostatics/scafacos_impl.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "grid.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

#include <cassert>
#include <memory>
#include <string>

std::shared_ptr<DipolarScafacos>
make_dipolar_scafacos(std::string const &method,
                      std::string const &parameters) {
  return std::make_shared<DipolarScafacosImpl>(comm_cart, method, parameters);
}

void DipolarScafacosImpl::update_particle_data() {
  positions.clear();
  dipoles.clear();

  for (auto const &p : cell_structure.local_particles()) {
    auto const pos = folded_position(p.pos(), box_geo);
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    auto const dip = p.calc_dip();
    dipoles.push_back(dip[0]);
    dipoles.push_back(dip[1]);
    dipoles.push_back(dip[2]);
  }
}

void DipolarScafacosImpl::update_particle_forces() const {
  if (positions.empty())
    return;

  auto it_potentials = potentials.begin();
  auto it_f = std::size_t{0ul};
  for (auto &p : cell_structure.local_particles()) {
    // The scafacos term "potential" here in fact refers to the magnetic
    // field. So, the torques are given by m \times B
    auto const dip = p.calc_dip();
    auto const t = vector_product(
        dip, Utils::Vector3d(Utils::Span<const double>(&*it_potentials, 3)));
    // The force is given by G m, where G is a matrix
    // which comes from the "fields" output of scafacos like this
    // 0 1 2
    // 1 3 4
    // 2 4 5
    // where the numbers refer to indices in the "field" output from scafacos
    auto const G = Utils::Matrix<double, 3, 3>{
        {fields[it_f + 0ul], fields[it_f + 1ul], fields[it_f + 2ul]},
        {fields[it_f + 1ul], fields[it_f + 3ul], fields[it_f + 4ul]},
        {fields[it_f + 2ul], fields[it_f + 4ul], fields[it_f + 5ul]}};
    auto const f = G * dip;

    // Add to particles
    p.force() += prefactor * f;
    p.torque() += prefactor * t;
    it_f += 6ul;
    it_potentials += 3;
  }

  /* Check that the particle number did not change */
  assert(it_potentials == potentials.end());
}

#endif // SCAFACOS_DIPOLES
