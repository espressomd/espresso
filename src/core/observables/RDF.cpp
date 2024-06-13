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
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/for_each_pair.hpp>
#include <utils/math/int_pow.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/range/combine.hpp>

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace Observables {
std::vector<double>
RDF::operator()(boost::mpi::communicator const &comm) const {
  auto const &local_particles_1 = fetch_particles(ids1());

  if (ids2().empty()) {
    return this->evaluate(comm, local_particles_1, {}, {});
  }

  auto const &local_particles_2 = fetch_particles(ids2());

  return this->evaluate(comm, local_particles_1, local_particles_2, {});
}

std::vector<double>
RDF::evaluate(boost::mpi::communicator const &comm,
              ParticleReferenceRange const &local_particles_1,
              ParticleReferenceRange const &local_particles_2,
              const ParticleObservables::traits<Particle> &traits) const {
  auto const positions_1 = detail::get_all_particle_positions(
      comm, local_particles_1, ids1(), traits, false);
  auto const positions_2 = detail::get_all_particle_positions(
      comm, local_particles_2, ids2(), traits, false);

  if (comm.rank() != 0) {
    return {};
  }

  auto const &box_geo = *System::get_system().box_geo;
  auto const bin_width = (max_r - min_r) / static_cast<double>(n_r_bins);
  auto const inv_bin_width = 1.0 / bin_width;
  std::vector<double> res(n_r_bins, 0.0);
  long int cnt = 0;
  auto op = [this, inv_bin_width, &cnt, &res, &box_geo](auto const &pos1,
                                                        auto const &pos2) {
    auto const dist = box_geo.get_mi_vector(pos1, pos2).norm();
    if (dist > min_r && dist < max_r) {
      auto const ind =
          static_cast<int>(std::floor((dist - min_r) * inv_bin_width));
      res[ind]++;
    }
    cnt++;
  };

  if (local_particles_2.empty()) {
    Utils::for_each_pair(positions_1, op);
  } else {
    auto const combine_1 = boost::combine(ids1(), positions_1);
    auto const combine_2 = boost::combine(ids2(), positions_2);

    auto op2 = [&op](auto const &it1, auto const &it2) {
      auto const &[id1, pos1] = it1;
      auto const &[id2, pos2] = it2;

      op(pos1, pos2);
    };

    auto cmp = [](auto const &it1, auto const &it2) {
      auto const &[id1, pos1] = it1;
      auto const &[id2, pos2] = it2;

      return id1 != id2;
    };
    Utils::for_each_cartesian_pair_if(combine_1, combine_2, op2, cmp);
  }
  if (cnt == 0)
    return res;
  // normalization
  auto const volume = box_geo.volume();
  for (std::size_t i = 0u; i < res.size(); ++i) {
    auto const r_in = static_cast<double>(i) * bin_width + min_r;
    auto const r_out = r_in + bin_width;
    auto const bin_volume =
        (4. / 3.) * std::numbers::pi *
        (Utils::int_pow<3>(r_out) - Utils::int_pow<3>(r_in));
    res[i] *= volume / (bin_volume * static_cast<double>(cnt));
  }

  return res;
}
} // namespace Observables
