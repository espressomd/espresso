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
/** \file
 *  Statistical tools to analyze simulations.
 *
 *  The corresponding header file is statistics.hpp.
 */

#include "analysis/statistics.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace boost::serialization {
/** @brief Serialize @c std::tuple. */
template <typename Archive, typename... T>
void serialize(Archive &ar, std::tuple<T...> &pack, unsigned int const) {
  std::apply([&](auto &...item) { ((ar & item), ...); }, pack);
}
} // namespace boost::serialization

/** @brief Decay a tuple of only 1 type to that type. */
template <typename... F> struct DecayTupleResult {
  using type = std::tuple<std::invoke_result_t<F, Particle const &>...>;
};

template <typename F> struct DecayTupleResult<F> {
  using type = std::invoke_result_t<F, Particle const &>;
};

/**
 * @brief Gather particle traits to MPI rank 0.
 * When only one trait is requested, return a vector of type T.
 * When multiple traits are requested, return a vector of tuples.
 */
template <class... Trait>
static auto gather_traits_for_types(System::System const &system,
                                    std::vector<int> const &p_types,
                                    Trait &&...trait) {
  std::vector<typename DecayTupleResult<Trait...>::type> buffer{};

  for (auto const &p : system.cell_structure->local_particles()) {
    if (Utils::contains(p_types, p.type())) {
      buffer.emplace_back(trait(p)...);
    }
  }

  Utils::Mpi::gather_buffer(buffer, ::comm_cart, 0);
  if (::comm_cart.rank() != 0) {
    buffer.clear();
  }
  return buffer;
}

double mindist(System::System const &system, std::vector<int> const &set1,
               std::vector<int> const &set2) {
  using Utils::contains;

  std::vector<int> buf_type{};
  std::vector<Utils::Vector3d> buf_pos{};

  auto const &box_geo = *system.box_geo;
  auto const accept_all = set1.empty() or set2.empty();
  for (auto const &p : system.cell_structure->local_particles()) {
    if (accept_all or contains(set1, p.type()) or contains(set2, p.type())) {
      buf_type.emplace_back(p.type());
      buf_pos.emplace_back(box_geo.unfolded_position(p.pos(), p.image_box()));
    }
  }

  Utils::Mpi::gather_buffer(buf_type, ::comm_cart, 0);
  Utils::Mpi::gather_buffer(buf_pos, ::comm_cart, 0);
  if (::comm_cart.rank() != 0) {
    buf_type.clear();
    buf_pos.clear();
  }

  auto mindist_sq = std::numeric_limits<double>::infinity();
  for (std::size_t j = 0ul; j < buf_type.size(); ++j) {
    /* check which sets particle j belongs to (bit 0: set1, bit1: set2) */
    auto in_set = 0;
    if (set1.empty() || contains(set1, buf_type[j]))
      in_set = 1;
    if (set2.empty() || contains(set2, buf_type[j]))
      in_set |= 2;
    assert(in_set != 0);

    for (auto i = j + 1ul; i != buf_type.size(); ++i) {
      /* accept a pair if particle j is in set1 and particle i in set2 or vice
       * versa. */
      if (((in_set & 1) and (set2.empty() or contains(set2, buf_type[i]))) or
          ((in_set & 2) and (set1.empty() or contains(set1, buf_type[i])))) {
        mindist_sq = std::min(
            mindist_sq, box_geo.get_mi_vector(buf_pos[j], buf_pos[i]).norm2());
      }
    }
  }

  return std::sqrt(mindist_sq);
}

Utils::Vector3d calc_linear_momentum(System::System const &system,
                                     bool include_particles,
                                     bool include_lbfluid) {
  Utils::Vector3d momentum{};
  if (include_particles) {
    auto const particles = system.cell_structure->local_particles();
    momentum =
        std::accumulate(particles.begin(), particles.end(), Utils::Vector3d{},
                        [](Utils::Vector3d const &m, Particle const &p) {
                          return m + p.mass() * p.v();
                        });
  }
  if (include_lbfluid and system.lb.is_solver_set()) {
    momentum += system.lb.get_momentum() * system.lb.get_lattice_speed();
  }
  return momentum;
}

Utils::Vector3d center_of_mass(System::System const &system, int p_type) {
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;
  Utils::Vector3d local_com{};
  double local_mass = 0.;

  for (auto const &p : cell_structure.local_particles()) {
    if ((p.type() == p_type or p_type == -1) and not p.is_virtual()) {
      local_com += box_geo.unfolded_position(p.pos(), p.image_box()) * p.mass();
      local_mass += p.mass();
    }
  }
  Utils::Vector3d com{};
  double mass = 1.; // placeholder value to avoid division by zero
  boost::mpi::reduce(::comm_cart, local_com, com, std::plus<>(), 0);
  boost::mpi::reduce(::comm_cart, local_mass, mass, std::plus<>(), 0);
  return com / mass;
}

Utils::Vector3d angular_momentum(System::System const &system, int p_type) {
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;
  Utils::Vector3d am{};

  for (auto const &p : cell_structure.local_particles()) {
    if ((p.type() == p_type or p_type == -1) and not p.is_virtual()) {
      auto const pos = box_geo.unfolded_position(p.pos(), p.image_box());
      am += p.mass() * vector_product(pos, p.v());
    }
  }
  return am;
}

Utils::Vector9d gyration_tensor(System::System const &system,
                                std::vector<int> const &p_types) {
  auto const &box_geo = *system.box_geo;
  auto const trait_pos = [&box_geo](Particle const &p) {
    return box_geo.unfolded_position(p.pos(), p.image_box());
  };
  auto const buf_pos = gather_traits_for_types(system, p_types, trait_pos);

  Utils::Vector9d mat{};
  if (::comm_cart.rank() == 0) {
    auto const center =
        std::accumulate(buf_pos.begin(), buf_pos.end(), Utils::Vector3d{}) /
        static_cast<double>(buf_pos.size());
    // compute covariance matrix
    for (unsigned int i = 0u; i < 3u; ++i) {
      for (unsigned int j = 0u; j < 3u; ++j) {
        if (i > j) {
          mat[i * 3u + j] = mat[j * 3u + i];
        } else {
          mat[i * 3u + j] = std::accumulate(
              buf_pos.begin(), buf_pos.end(), 0.,
              [i, j, &center](double acc, Utils::Vector3d const &pos) {
                return acc + (pos[i] - center[i]) * (pos[j] - center[j]);
              });
        }
      }
    }
    mat /= static_cast<double>(buf_pos.size());
  }
  return mat;
}

Utils::Vector9d moment_of_inertia_matrix(System::System const &system,
                                         int p_type) {
  auto const &box_geo = *system.box_geo;
  auto const &cell_structure = *system.cell_structure;
  Utils::Vector9d mat{};
  auto com = center_of_mass(system, p_type);
  boost::mpi::broadcast(::comm_cart, com, 0);

  for (auto const &p : cell_structure.local_particles()) {
    if (p.type() == p_type and not p.is_virtual()) {
      auto const pos = box_geo.unfolded_position(p.pos(), p.image_box()) - com;
      auto const mass = p.mass();
      mat[0] += mass * (pos[1] * pos[1] + pos[2] * pos[2]);
      mat[4] += mass * (pos[0] * pos[0] + pos[2] * pos[2]);
      mat[8] += mass * (pos[0] * pos[0] + pos[1] * pos[1]);
      mat[1] -= mass * (pos[0] * pos[1]);
      mat[2] -= mass * (pos[0] * pos[2]);
      mat[5] -= mass * (pos[1] * pos[2]);
    }
  }
  /* use symmetry */
  mat[3] = mat[1];
  mat[6] = mat[2];
  mat[7] = mat[5];
  return mat;
}

std::vector<int> nbhood(System::System const &system,
                        Utils::Vector3d const &pos, double dist) {
  std::vector<int> buf_pid{};
  auto const dist_sq = dist * dist;
  auto const &box_geo = *system.box_geo;

  for (auto const &p : system.cell_structure->local_particles()) {
    auto const r_sq = box_geo.get_mi_vector(pos, p.pos()).norm2();
    if (r_sq < dist_sq) {
      buf_pid.push_back(p.id());
    }
  }

  Utils::Mpi::gather_buffer(buf_pid, ::comm_cart, 0);
  if (::comm_cart.rank() != 0) {
    buf_pid.clear();
  }

  return buf_pid;
}

std::vector<std::vector<double>>
calc_part_distribution(System::System const &system,
                       std::vector<int> const &p1_types,
                       std::vector<int> const &p2_types, double r_min,
                       double r_max, int r_bins, bool log_flag, bool int_flag) {

  auto const &box_geo = *system.box_geo;
  auto const trait_id = [](Particle const &p) { return p.id(); };
  auto const trait_pos = [&box_geo](Particle const &p) {
    return box_geo.unfolded_position(p.pos(), p.image_box());
  };
  auto const buf1 =
      gather_traits_for_types(system, p1_types, trait_id, trait_pos);
  auto const buf2 =
      gather_traits_for_types(system, p2_types, trait_id, trait_pos);
  auto const r_max2 = Utils::sqr(r_max);
  auto const r_min2 = Utils::sqr(r_min);
  auto const start_dist2 = Utils::sqr(r_max + 1.);
  auto const inv_bin_width =
      (log_flag) ? static_cast<double>(r_bins) / std::log(r_max / r_min)
                 : static_cast<double>(r_bins) / (r_max - r_min);

  long cnt = 0;
  double low = 0.0;
  std::vector<double> distribution(r_bins);

  for (auto const &[pid1, pos1] : buf1) {
    auto min_dist2 = start_dist2;
    /* particle loop: p2_types */
    for (auto const &[pid2, pos2] : buf2) {
      if (pid1 != pid2) {
        auto const act_dist2 = box_geo.get_mi_vector(pos1, pos2).norm2();
        if (act_dist2 < min_dist2) {
          min_dist2 = act_dist2;
        }
      }
    }
    if (min_dist2 <= r_max2) {
      if (min_dist2 >= r_min2) {
        auto const min_dist = std::sqrt(min_dist2);
        /* calculate bin index */
        auto const ind = static_cast<int>(
            ((log_flag) ? std::log(min_dist / r_min) : (min_dist - r_min)) *
            inv_bin_width);
        if (ind >= 0 and ind < r_bins) {
          distribution[ind] += 1.0;
        }
      } else {
        low += 1.0;
      }
    }
    cnt++;
  }

  if (cnt != 0) {
    // normalization
    low /= static_cast<double>(cnt);
    for (int i = 0; i < r_bins; i++) {
      distribution[i] /= static_cast<double>(cnt);
    }

    // integration
    if (int_flag) {
      distribution[0] += low;
      for (int i = 0; i < r_bins - 1; i++)
        distribution[i + 1] += distribution[i];
    }
  }

  std::vector<double> radii(r_bins);
  if (log_flag) {
    auto const log_fac = std::pow(r_max / r_min, 1. / r_bins);
    radii[0] = r_min * std::sqrt(log_fac);
    for (int i = 1; i < r_bins; ++i) {
      radii[i] = radii[i - 1] * log_fac;
    }
  } else {
    auto const bin_width = (r_max - r_min) / static_cast<double>(r_bins);
    for (int i = 0; i < r_bins; ++i) {
      radii[i] = r_min + bin_width / 2. + static_cast<double>(i) * bin_width;
    }
  }

  return {radii, distribution};
}

std::vector<std::vector<double>>
structure_factor(System::System const &system, std::vector<int> const &p_types,
                 int order) {
  auto const &box_geo = *system.box_geo;
  auto const trait_pos = [&box_geo](Particle const &p) {
    return box_geo.unfolded_position(p.pos(), p.image_box());
  };
  auto const buf_pos = gather_traits_for_types(system, p_types, trait_pos);
  auto const order_sq = Utils::sqr(static_cast<std::size_t>(order));
  auto const twoPI_L = 2. * std::numbers::pi * system.box_geo->length_inv()[0];
  std::vector<double> ff(2ul * order_sq + 1ul);
  std::vector<double> wavevectors;
  std::vector<double> intensities;

  if (::comm_cart.rank() == 0) {
    for (int i = 0; i <= order; i++) {
      for (int j = -order; j <= order; j++) {
        for (int k = -order; k <= order; k++) {
          auto const n = i * i + j * j + k * k;
          if ((static_cast<std::size_t>(n) <= order_sq) && (n >= 1)) {
            double C_sum = 0.0, S_sum = 0.0;
            for (auto const &pos : buf_pos) {
              auto const qr = twoPI_L * (Utils::Vector3i{{i, j, k}} * pos);
              C_sum += cos(qr);
              S_sum += sin(qr);
            }
            ff[2 * n - 2] += C_sum * C_sum + S_sum * S_sum;
            ff[2 * n - 1]++;
          }
        }
      }
    }

    std::size_t length = 0;
    for (std::size_t qi = 0; qi < order_sq; qi++) {
      if (ff[2 * qi + 1] != 0) {
        ff[2 * qi] /= static_cast<double>(buf_pos.size()) * ff[2 * qi + 1];
        ++length;
      }
    }

    wavevectors.resize(length);
    intensities.resize(length);

    int cnt = 0;
    for (std::size_t i = 0; i < order_sq; i++) {
      if (ff[2 * i + 1] != 0) {
        wavevectors[cnt] = twoPI_L * std::sqrt(static_cast<double>(i + 1));
        intensities[cnt] = ff[2 * i];
        cnt++;
      }
    }
  }

  return {std::move(wavevectors), std::move(intensities)};
}
