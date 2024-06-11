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
 *  This file contains everything needed to create a start-up configuration
 *  of (possibly charged) polymer chains with counterions and salt molecules,
 *  assigning velocities to the particles and cross-linking the polymers if
 *  necessary.
 *
 *  The corresponding header file is polymer.hpp.
 */

#include "polymer.hpp"

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/Constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"
#include "random.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/vec_rotate.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

template <class RNG>
static Utils::Vector3d random_position(BoxGeometry const &box_geo, RNG &rng) {
  Utils::Vector3d v;
  for (auto i = 0u; i < 3u; ++i)
    v[i] = box_geo.length()[i] * rng();
  return v;
}

template <class RNG> static Utils::Vector3d random_unit_vector(RNG &rng) {
  Utils::Vector3d v;
  double const phi = acos(1. - 2. * rng());
  double const theta = 2. * Utils::pi() * rng();
  v[0] = sin(phi) * cos(theta);
  v[1] = sin(phi) * sin(theta);
  v[2] = cos(phi);
  v /= v.norm();
  return v;
}

/** Determines whether a given position @p pos is valid, i.e., it doesn't
 *  collide with existing or buffered particles, nor with existing constraints
 *  (if @c respect_constraints).
 *  @param pos                   the trial position in question
 *  @param positions             buffered positions to respect
 *  @param cell_structure        existing particles to respect
 *  @param box_geo               Box geometry
 *  @param min_distance          threshold for the minimum distance between
 *                               trial position and buffered/existing particles
 *  @param respect_constraints   whether to respect constraints
 *  @return true if valid position, false if not.
 */
static bool
is_valid_position(Utils::Vector3d const &pos,
                  std::vector<std::vector<Utils::Vector3d>> const &positions,
                  CellStructure const &cell_structure,
                  BoxGeometry const &box_geo, double const min_distance,
                  int const respect_constraints) {

  struct reduce_min {
    auto operator()(double const a, double const b) const {
      return std::min(a, b);
    }
  };

  // check if constraint is violated
  if (respect_constraints) {
    Utils::Vector3d const folded_pos = box_geo.folded_position(pos);

    for (auto &c : Constraints::constraints) {
      auto cs =
          std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
      if (cs) {
        double d;
        Utils::Vector3d v;

        cs->calc_dist(folded_pos, d, v);

        if (d <= 0) {
          return false;
        }
      }
    }
  }

  if (min_distance > 0.) {
    // check for collision with existing particles
    auto local_mindist_sq = std::numeric_limits<double>::infinity();
    for (auto const &p : cell_structure.local_particles()) {
      auto const d = box_geo.get_mi_vector(pos, p.pos());
      local_mindist_sq = std::min(local_mindist_sq, d.norm2());
    }
    auto const global_mindist_sq =
        boost::mpi::all_reduce(::comm_cart, local_mindist_sq, reduce_min{});
    if (std::sqrt(global_mindist_sq) < min_distance) {
      return false;
    }

    auto const min_distance_sq = min_distance * min_distance;
    for (auto const &p : positions) {
      for (auto const &m : p) {
        if (box_geo.get_mi_vector(pos, m).norm2() < min_distance_sq) {
          return false;
        }
      }
    }
  }
  return true;
}

std::vector<std::vector<Utils::Vector3d>>
draw_polymer_positions(System::System const &system, int const n_polymers,
                       int const beads_per_chain, double const bond_length,
                       std::vector<Utils::Vector3d> const &start_positions,
                       double const min_distance, int const max_tries,
                       int const use_bond_angle, double const bond_angle,
                       int const respect_constraints, int const seed) {

  auto const &box_geo = *system.box_geo;
  auto rng = [mt = Random::mt19937(static_cast<unsigned>(seed)),
              dist = std::uniform_real_distribution<double>(
                  0.0, 1.0)]() mutable { return dist(mt); };

  std::vector<std::vector<Utils::Vector3d>> positions(n_polymers);
  for (auto &p : positions) {
    p.reserve(beads_per_chain);
  }

  auto is_valid_pos = [&](Utils::Vector3d const &pos) {
    return is_valid_position(pos, positions, *system.cell_structure, box_geo,
                             min_distance, respect_constraints);
  };

  for (std::size_t p = 0; p < start_positions.size(); p++) {
    if (is_valid_pos(start_positions[p])) {
      positions[p].push_back(start_positions[p]);
    } else {
      throw std::runtime_error("Invalid start positions.");
    }
  }

  /* Draw a monomer position, obeying angle and starting position
   * constraints where appropriate. */
  auto draw_monomer_position = [&](int p, int m) {
    if (m == 0) {
      return (p < start_positions.size()) ? start_positions[p]
                                          : random_position(box_geo, rng);
    }

    if (not use_bond_angle or m < 2) {
      return positions[p][m - 1] + bond_length * random_unit_vector(rng);
    }

    auto const last_vec = positions[p][m - 1] - positions[p][m - 2];
    return positions[p][m - 1] +
           Utils::vec_rotate(vector_product(last_vec, random_unit_vector(rng)),
                             bond_angle, -last_vec);
  };

  /* Try up to max_tries times to draw a valid position */
  auto draw_valid_monomer_position =
      [&](int p, int m) -> std::optional<Utils::Vector3d> {
    for (auto i = 0; i < max_tries; i++) {
      auto const trial_pos = draw_monomer_position(p, m);
      if (is_valid_pos(trial_pos)) {
        return trial_pos;
      }
    }

    return std::nullopt;
  };

  // create remaining monomers' positions by backtracking.
  for (int p = 0; p < n_polymers; ++p) {
    for (int attempts_poly = 0; attempts_poly < max_tries; attempts_poly++) {
      int rejections = 0;
      while (positions[p].size() < beads_per_chain) {
        auto pos = draw_valid_monomer_position(
            p, static_cast<int>(positions[p].size()));

        if (pos) {
          /* Move on one position */
          positions[p].push_back(*pos);
        } else if (not positions[p].empty()) {
          /* Go back one position and try again */
          positions[p].pop_back();
          rejections++;
          if (rejections > max_tries) {
            /* Give up for this try. */
            break;
          }
        } else {
          /* Give up for this try. */
          break;
        }
      }

      /* If the polymer has not full length, we  try again. Otherwise we can
       * move on to the next polymer. */
      if (positions[p].size() == beads_per_chain) {
        break;
      }
    }

    /* We did not get a complete polymer, but have exceeded the maximal
     * number of tries, which means failure. */
    if (positions[p].size() < beads_per_chain)
      throw std::runtime_error("Failed to create polymer positions.");
  }
  return positions;
}
