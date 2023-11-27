/*
 * Copyright (C) 2019-2022 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLEDIHEDRALS_HPP
#define OBSERVABLES_PARTICLEDIHEDRALS_HPP

#include "BoxGeometry.hpp"
#include "PidObservable.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** Calculate dihedral angles between particles in a polymer.
 *  For @f$ n @f$ bonded particles, return the @f$ n-3 @f$ dihedral angles
 *  along the chain, in radians.
 *
 *  The sign of the dihedral angles follow the IUPAC nomenclature for the
 *  Newman projection, specifically section "Torsion Angle" pages 2220-2221
 *  in @cite moss96a.
 *
 */
class BondDihedrals : public PidObservable {
public:
  using PidObservable::PidObservable;
  explicit BondDihedrals(std::vector<int> ids) : PidObservable(std::move(ids)) {
    if (this->ids().size() < 4)
      throw std::runtime_error("At least 4 particles are required");
  }

  std::vector<double>
  evaluate(boost::mpi::communicator const &comm,
           ParticleReferenceRange const &local_particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    auto const positions_sorted = detail::get_all_particle_positions(
        comm, local_particles, ids(), traits, false);

    if (comm.rank() != 0) {
      return {};
    }

    auto const &box_geo = *System::get_system().box_geo;
    std::vector<double> res(n_values());
    auto v1 = box_geo.get_mi_vector(positions_sorted[1], positions_sorted[0]);
    auto v2 = box_geo.get_mi_vector(positions_sorted[2], positions_sorted[1]);
    auto c1 = Utils::vector_product(v1, v2);
    for (std::size_t i = 0, end = n_values(); i < end; i++) {
      auto v3 = box_geo.get_mi_vector(positions_sorted[i + 3],
                                      positions_sorted[i + 2]);
      auto c2 = vector_product(v2, v3);
      /* the 2-argument arctangent returns an angle in the range [-pi, pi] that
       * allows for an unambiguous determination of the 4th particle position */
      res[i] = atan2((Utils::vector_product(c1, c2) * v2) / v2.norm(), c1 * c2);
      v1 = v2;
      v2 = v3;
      c1 = c2;
    }
    return res;
  }
  std::vector<std::size_t> shape() const override {
    assert(ids().size() >= 3);
    return {ids().size() - 3};
  }
};

} // Namespace Observables

#endif
