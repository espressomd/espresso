/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef CORE_SHORT_RANGE_HPP
#define CORE_SHORT_RANGE_HPP

#include "algorithm/for_each_pair.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <boost/iterator/indirect_iterator.hpp>
#include <profiler/profiler.hpp>

#include <utility>

/**
 * @brief Distance vector and length handed to pair kernels.
 */
struct Distance {
  explicit Distance(Utils::Vector3d const &vec21)
      : vec21(vec21), dist2(vec21.norm2()) {}

  Utils::Vector3d vec21;
  const double dist2;
};

namespace detail {
struct MinimalImageDistance {
  const BoxGeometry box;

  Distance operator()(Particle const &p1, Particle const &p2) const {
    return Distance(get_mi_vector(p1.r.p, p2.r.p, box));
  }
};

} // namespace detail

template <typename ParticleKernel, typename PairKernel>
void short_range_loop(ParticleKernel &&particle_kernel,
                      PairKernel &&pair_kernel) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  assert(get_resort_particles() =Cells::RESORT_NONE);

  auto first = boost::make_indirect_iterator(local_cells.begin());
  auto last = boost::make_indirect_iterator(local_cells.end());

#ifdef ELECTROSTATICS
  auto const coulomb_cutoff = Coulomb::cutoff(box_geo.length());
#else
  auto const coulomb_cutoff = INACTIVE_CUTOFF;
#endif

#ifdef DIPOLES
  auto const dipole_cutoff = Dipole::cutoff(box_geo.length());
#else
  auto const dipole_cutoff = INACTIVE_CUTOFF;
#endif

  ParticleKernel &&particleKernel =
      std::forward<ParticleKernel>(particle_kernel);
  PairKernel &&pairKernel = std::forward<PairKernel>(pair_kernel);
  VerletCriterion &&verletCriterion =
      VerletCriterion{skin, max_cut, coulomb_cutoff, dipole_cutoff,
                      collision_detection_cutoff()};

  Algorithm::for_each_pair(
      first, last, std::forward<ParticleKernel>(particleKernel),
      std::forward<PairKernel>(pairKernel), detail::MinimalImageDistance{box_geo},
      std::forward<VerletCriterion>(verletCriterion),
      cell_structure.use_verlet_list, rebuild_verletlist);

  rebuild_verletlist = 0;
}

#endif
