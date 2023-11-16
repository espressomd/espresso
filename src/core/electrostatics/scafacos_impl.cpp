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

#ifdef SCAFACOS

#include "electrostatics/scafacos.hpp"
#include "electrostatics/scafacos_impl.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"
#include "tuning.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <boost/range/algorithm/min_element.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <string>

std::shared_ptr<CoulombScafacos>
make_coulomb_scafacos(std::string const &method,
                      std::string const &parameters) {
  return std::make_shared<CoulombScafacosImpl>(comm_cart, method, parameters);
}

void CoulombScafacosImpl::update_particle_data() {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;

  positions.clear();
  charges.clear();

  for (auto const &p : cell_structure.local_particles()) {
    auto const pos = box_geo.folded_position(p.pos());
    positions.push_back(pos[0]);
    positions.push_back(pos[1]);
    positions.push_back(pos[2]);
    charges.push_back(p.q());
  }
}

void CoulombScafacosImpl::update_particle_forces() const {
  if (positions.empty())
    return;

  auto const &cell_structure = *get_system().cell_structure;

  auto it_fields = fields.begin();
  for (auto &p : cell_structure.local_particles()) {
    p.force() += prefactor * p.q() *
                 Utils::Vector3d(Utils::Span<const double>(&*it_fields, 3));
    it_fields += 3;
  }

  /* Check that the particle number did not change */
  assert(it_fields == fields.end());
}

double CoulombScafacosImpl::time_r_cut(double r_cut) {
  set_r_cut_and_tune(r_cut);
  auto &system = get_system();
  return benchmark_integration_step(system, 10);
}

void CoulombScafacosImpl::tune_r_cut() {
  auto constexpr convergence_threshold = 1e-3;
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const skin = system.get_verlet_skin();

  auto const min_box_l = *boost::min_element(box_geo.length());
  auto const min_local_box_l = *boost::min_element(local_geo.length());

  /* The bisection code breaks down when r_min < 1 for several methods
   * (e.g. p2nfft, p3m, ewald) if the mesh size is not fixed (ScaFaCoS
   * either hangs or allocates too much memory) */
  auto r_min = 1.0;
  auto r_max = std::min(min_local_box_l, min_box_l / 2.0) - skin;
  assert(r_max >= r_min);
  auto t_min = 0.0;
  auto t_max = std::numeric_limits<double>::max();
  auto r_opt = -1.0;

  /* Run bisection */
  while (std::fabs(r_min - r_max) > convergence_threshold) {
    r_opt = (r_max + r_min) / 2.;
    auto const dr = 0.5 * (r_max - r_min);
    auto const t_mid = time_r_cut(r_min + dr);
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
  assert(r_opt >= 0.);
  set_r_cut(r_opt);
}

void CoulombScafacosImpl::tune_impl() {
  update_particle_data();

  // Check whether we have to do a bisection for the short-range cutoff
  // Check if there is a user-supplied cutoff
  if (ScafacosContext::get_near_field_delegation() and
      ScafacosContext::r_cut() <= 0.0) {
    tune_r_cut();
  } else {
    // ESPResSo is not affected by a short-range cutoff -> tune in parallel
    ScafacosContext::tune(charges, positions);
  }
  get_system().on_coulomb_change();
}

#endif // SCAFACOS
