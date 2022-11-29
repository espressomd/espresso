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
#ifndef OBSERVABLES_BONDANGLES_HPP
#define OBSERVABLES_BONDANGLES_HPP

#include "BoxGeometry.hpp"
#include "PidObservable.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** Calculate bond angles between particles in a polymer.
 *  For @f$ n @f$ bonded particles, return the @f$ n-2 @f$ angles along the
 *  chain, in radians.
 */
class BondAngles : public PidObservable {
public:
  using PidObservable::PidObservable;
  explicit BondAngles(std::vector<int> ids) : PidObservable(std::move(ids)) {
    if (this->ids().size() < 3)
      throw std::runtime_error("At least 3 particles are required");
  }

  std::vector<double>
  evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    std::vector<double> res(n_values());
    auto v1 = box_geo.get_mi_vector(traits.position(particles[1]),
                                    traits.position(particles[0]));
    auto n1 = v1.norm();
    for (std::size_t i = 0, end = n_values(); i < end; i++) {
      auto v2 = box_geo.get_mi_vector(traits.position(particles[i + 2]),
                                      traits.position(particles[i + 1]));
      auto const n2 = v2.norm();
      auto const cosine =
          std::clamp((v1 * v2) / (n1 * n2), -TINY_COS_VALUE, TINY_COS_VALUE);
      /* to reduce computational time, after calculating an angle ijk, the
       * vector r_ij takes the value r_jk, but to orient it correctly, it has
       * to be multiplied -1; it's cheaper to do this operation on a double
       * than on a vector of doubles
       */
      res[i] = acos(-cosine);
      v1 = v2;
      n1 = n2;
    }
    return res;
  }
  std::vector<std::size_t> shape() const override {
    assert(ids().size() >= 2);
    return {ids().size() - 2};
  }
};

} // Namespace Observables

#endif
