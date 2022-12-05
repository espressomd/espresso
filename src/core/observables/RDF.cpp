/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include "RDF.hpp"

#include "BoxGeometry.hpp"
#include "fetch_particles.hpp"
#include "grid.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/for_each_pair.hpp>
#include <utils/math/int_pow.hpp>

#include <boost/range/algorithm/transform.hpp>

#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace Observables {
std::vector<double> RDF::operator()() const {
  std::vector<Particle> particles1 = fetch_particles(ids1());
  std::vector<const Particle *> particles_ptrs1(particles1.size());
  boost::transform(particles1, particles_ptrs1.begin(),
                   [](auto const &p) { return std::addressof(p); });

  if (ids2().empty()) {
    return this->evaluate(particles_ptrs1, {});
  }

  std::vector<Particle> particles2 = fetch_particles(ids2());
  std::vector<const Particle *> particles_ptrs2(particles2.size());
  boost::transform(particles2, particles_ptrs2.begin(),
                   [](auto const &p) { return std::addressof(p); });
  return this->evaluate(particles_ptrs1, particles_ptrs2);
}

std::vector<double>
RDF::evaluate(Utils::Span<const Particle *const> particles1,
              Utils::Span<const Particle *const> particles2) const {
  auto const bin_width = (max_r - min_r) / static_cast<double>(n_r_bins);
  auto const inv_bin_width = 1.0 / bin_width;
  std::vector<double> res(n_values(), 0.0);
  long int cnt = 0;
  auto op = [this, inv_bin_width, &cnt, &res](const Particle *const p1,
                                              const Particle *const p2) {
    auto const dist = box_geo.get_mi_vector(p1->pos(), p2->pos()).norm();
    if (dist > min_r && dist < max_r) {
      auto const ind =
          static_cast<int>(std::floor((dist - min_r) * inv_bin_width));
      res[ind]++;
    }
    cnt++;
  };
  if (particles2.empty()) {
    Utils::for_each_pair(particles1, op);
  } else {
    auto cmp = std::not_equal_to<const Particle *const>();
    Utils::for_each_cartesian_pair_if(particles1, particles2, op, cmp);
  }
  if (cnt == 0)
    return res;
  // normalization
  auto const volume = box_geo.volume();
  for (int i = 0; i < n_r_bins; ++i) {
    auto const r_in = i * bin_width + min_r;
    auto const r_out = r_in + bin_width;
    auto const bin_volume =
        (4.0 / 3.0) * Utils::pi() *
        (Utils::int_pow<3>(r_out) - Utils::int_pow<3>(r_in));
    res[i] *= volume / (bin_volume * static_cast<double>(cnt));
  }

  return res;
}
} // namespace Observables
