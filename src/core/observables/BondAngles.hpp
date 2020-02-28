/*
 * Copyright (C) 2019 The ESPResSo project
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

#include "PidObservable.hpp"
#include <utils/Vector.hpp>

#include <cmath>
#include <vector>

namespace Observables {

/** Calculate bond angles between particles in a polymer.
 *  For @f$ n @f$ bonded particles, return the @f$ n-2 @f$ angles along the
 *  chain, in radians.
 */
class BondAngles : public PidObservable {
public:
  using PidObservable::PidObservable;

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::vector<double> res(n_values());
    auto v1 = get_mi_vector(particles[1]->r.p, particles[0]->r.p, box_geo);
    auto n1 = v1.norm();
    for (size_t i = 0, end = n_values(); i < end; i++) {
      auto v2 =
          get_mi_vector(particles[i + 2]->r.p, particles[i + 1]->r.p, box_geo);
      auto n2 = v2.norm();
      auto cosine = (v1 * v2) / (n1 * n2);
      // sanitize cosine value
      if (cosine > TINY_COS_VALUE)
        cosine = TINY_COS_VALUE;
      else if (cosine < -TINY_COS_VALUE)
        cosine = -TINY_COS_VALUE;
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
  std::vector<size_t> shape() const override { return {ids().size() - 2}; }
};

} // Namespace Observables

#endif
