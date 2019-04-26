/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file
    This file contains everything needed to create a start-up configuration
    of (partially charged) polymer chains with counterions and salt molecules,
    assigning velocities to the particles and cross-linking the polymers if
   necessary.

    The corresponding header file is polymer.hpp.
*/

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "PartCfg.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"
#include "debug.hpp"
#include "global.hpp"
#include "integrate.hpp"
#include "polymer.hpp"
#include "random.hpp"
#include "utils.hpp"

#include "utils/vec_rotate.hpp"

Utils::Vector3d random_position(const std::function<double()> &generate_rn) {
  Utils::Vector3d v;
  for (int i = 0; i < 3; ++i)
    v[i] = box_l[i] * generate_rn();
  return v;
}

Utils::Vector3d random_unit_vector(const std::function<double()> &generate_rn) {
  Utils::Vector3d v;
  const double phi = acos(1. - 2. * generate_rn());
  const double theta = 2. * Utils::pi() * generate_rn();
  v[0] = sin(phi) * cos(theta);
  v[1] = sin(phi) * sin(theta);
  v[2] = cos(phi);
  v /= v.norm();
  return v;
}

double mindist(PartCfg &partCfg, const Utils::Vector3d pos) {
  if (partCfg.size() == 0) {
    return std::min(std::min(box_l[0], box_l[1]), box_l[2]);
  }

  auto const mindist = std::accumulate(
      partCfg.begin(), partCfg.end(), std::numeric_limits<double>::infinity(),
      [&pos](double mindist, Particle const &p) {
        return std::min(mindist, get_mi_vector(pos, p.r.p).norm2());
      });

  if (mindist < std::numeric_limits<double>::infinity())
    return std::sqrt(mindist);
  return -1.0;
}

bool is_valid_position(
    const Utils::Vector3d *pos,
    const std::vector<std::vector<Utils::Vector3d>> *positions,
    PartCfg &partCfg, const double min_distance,
    const int respect_constraints) {
  const Utils::Vector3d folded_pos = folded_position(*pos);
  // check if constraint is violated
  if (respect_constraints) {
    for (auto &c : Constraints::constraints) {
      auto cs =
          std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
      if (cs) {
        double d;
        double v[3];

        cs->calc_dist(folded_pos, &d, v);

        if (d <= 0) {
          return false;
        }
      }
    }
  }

  if (min_distance > 0) {
    // check for collision with existing particles
    if (mindist(partCfg, *pos) < min_distance) {
      return false;
    }
    // check for collision with buffered positions
    double buff_mindist = std::numeric_limits<double>::infinity();
    double h;
    for (const auto p : *positions) {
      for (auto m : p) {
        h = (folded_position(*pos) - folded_position(m)).norm2();
        buff_mindist = std::min(h, buff_mindist);
      }
    }
    if (std::sqrt(buff_mindist) < min_distance) {
      return false;
    }
  }
  return true;
}

std::vector<std::vector<Utils::Vector3d>>
draw_polymer_positions(PartCfg &partCfg, const int n_polymers,
                       const int beads_per_chain, const double bond_length,
                       const std::vector<Utils::Vector3d> &start_positions,
                       const double min_distance, const int max_tries,
                       const int use_bond_angle, const double bond_angle,
                       const int respect_constraints, const int seed) {
  std::vector<std::vector<Utils::Vector3d>> positions(
      n_polymers, std::vector<Utils::Vector3d>(beads_per_chain));

  std::mt19937 mt(seed);
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  Utils::Vector3d trial_pos;
  int attempts_mono, attempts_poly;

  // make sure that if given, all starting positions are valid
  if ((not start_positions.empty()) and
      std::any_of(start_positions.begin(), start_positions.end(),
                  [&positions, &partCfg, min_distance,
                   respect_constraints](Utils::Vector3d v) {
                    return not is_valid_position(&v, &positions, partCfg,
                                                 min_distance,
                                                 respect_constraints);
                  }))
    throw std::runtime_error("Invalid start positions.");
  // use (if none given, random) starting positions for every first monomer
  for (int p = 0; p < n_polymers; ++p) {
    attempts_mono = 0;
    // first monomer for all polymers
    if (start_positions.empty()) {
      do {
        trial_pos = random_position([&]() { return dist(mt); });
        attempts_mono++;
      } while ((not is_valid_position(&trial_pos, &positions, partCfg,
                                      min_distance, respect_constraints)) and
               (attempts_mono < max_tries));
      if (attempts_mono == max_tries) {
        throw std::runtime_error("Failed to create polymer start positions.");
      }
      positions[p][0] = trial_pos;
    } else {
      positions[p][0] = start_positions[p];
    }
  }

  // create remaining monomers' positions
  for (int p = 0; p < n_polymers; ++p) {
    attempts_poly = 0;
    for (int m = 1; m < beads_per_chain; ++m) {
      attempts_mono = 0;
      do {
        if (m == 0) {
          // m == 0 is only true after a failed attempt to position a polymer
          trial_pos = random_position([&]() { return dist(mt); });
        } else if (not use_bond_angle or m < 2) {
          // random step, also necessary if angle is set, placing the second
          // bead
          trial_pos =
              positions[p][m - 1] +
              bond_length * random_unit_vector([&]() { return dist(mt); });
        } else {
          // use prescribed angle
          Utils::Vector3d last_vec = positions[p][m - 1] - positions[p][m - 2];
          trial_pos = positions[p][m - 1] +
                      Utils::vec_rotate(
                          vector_product(last_vec, random_unit_vector([&]() {
                                           return dist(mt);
                                         })),
                          bond_angle, -last_vec);
        }
        attempts_mono++;
      } while ((not is_valid_position(&trial_pos, &positions, partCfg,
                                      min_distance, respect_constraints)) and
               (attempts_mono < max_tries));

      if (attempts_mono == max_tries) {
        if (attempts_poly < max_tries) {
          // if start positions have to be respected: fail ...
          if (start_positions.empty()) {
            throw std::runtime_error("Failed to create polymer positions with "
                                     "given start positions.");
          }
          // ... otherwise retry to position the whole polymer
          attempts_poly++;
          m = -1;
        } else {
          // ... but only if max_tries has not been exceeded.
          throw std::runtime_error("Failed to create polymer positions.");
        }
      } else {
        positions[p][m] = trial_pos;
      }
    }
  }
  return positions;
}
