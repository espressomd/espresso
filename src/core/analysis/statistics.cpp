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

#include "Particle.hpp"
#include "cells.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "partCfg_global.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/contains.hpp>
#include <utils/math/sqr.hpp>

#include <cassert>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>

double mindist(PartCfg &partCfg, std::vector<int> const &set1,
               std::vector<int> const &set2) {
  using Utils::contains;

  auto mindist_sq = std::numeric_limits<double>::infinity();

  for (auto jt = partCfg.begin(); jt != partCfg.end(); ++jt) {
    /* check which sets particle j belongs to (bit 0: set1, bit1: set2) */
    auto in_set = 0u;
    if (set1.empty() || contains(set1, jt->type()))
      in_set = 1u;
    if (set2.empty() || contains(set2, jt->type()))
      in_set |= 2u;
    if (in_set == 0)
      continue;

    for (auto it = std::next(jt); it != partCfg.end(); ++it)
      /* accept a pair if particle j is in set1 and particle i in set2 or vice
       * versa. */
      if (((in_set & 1u) && (set2.empty() || contains(set2, it->type()))) ||
          ((in_set & 2u) && (set1.empty() || contains(set1, it->type()))))
        mindist_sq = std::min(
            mindist_sq, box_geo.get_mi_vector(jt->pos(), it->pos()).norm2());
  }

  return std::sqrt(mindist_sq);
}

Utils::Vector3d calc_linear_momentum(bool include_particles,
                                     bool include_lbfluid) {
  Utils::Vector3d momentum{};
  if (include_particles) {
    auto const particles = cell_structure.local_particles();
    momentum =
        std::accumulate(particles.begin(), particles.end(), Utils::Vector3d{},
                        [](Utils::Vector3d const &m, Particle const &p) {
                          return m + p.mass() * p.v();
                        });
  }
  if (include_lbfluid and lattice_switch != ActiveLB::NONE) {
    momentum += LB::calc_fluid_momentum() * LB::get_lattice_speed();
  }
  return momentum;
}

Utils::Vector3d center_of_mass(PartCfg &partCfg, int p_type) {
  Utils::Vector3d com{};
  double mass = 0.0;

  for (auto const &p : partCfg) {
    if ((p.type() == p_type or p_type == -1) and not p.is_virtual()) {
      com += p.pos() * p.mass();
      mass += p.mass();
    }
  }
  com /= mass;
  return com;
}

Utils::Vector3d angular_momentum(PartCfg &partCfg, int p_type) {
  Utils::Vector3d am{};

  for (auto const &p : partCfg) {
    if ((p.type() == p_type or p_type == -1) and not p.is_virtual()) {
      am += p.mass() * vector_product(p.pos(), p.v());
    }
  }
  return am;
}

Utils::Vector9d moment_of_inertia_matrix(PartCfg &partCfg, int p_type) {
  Utils::Vector9d mat{};
  auto const com = center_of_mass(partCfg, p_type);
  for (auto const &p : partCfg) {
    if (p.type() == p_type and (not p.is_virtual())) {
      auto const p1 = p.pos() - com;
      auto const mass = p.mass();
      mat[0] += mass * (p1[1] * p1[1] + p1[2] * p1[2]);
      mat[4] += mass * (p1[0] * p1[0] + p1[2] * p1[2]);
      mat[8] += mass * (p1[0] * p1[0] + p1[1] * p1[1]);
      mat[1] -= mass * (p1[0] * p1[1]);
      mat[2] -= mass * (p1[0] * p1[2]);
      mat[5] -= mass * (p1[1] * p1[2]);
    }
  }
  /* use symmetry */
  mat[3] = mat[1];
  mat[6] = mat[2];
  mat[7] = mat[5];
  return mat;
}

std::vector<int> nbhood(PartCfg &partCfg, Utils::Vector3d const &pos,
                        double dist) {
  std::vector<int> ids;
  auto const dist_sq = dist * dist;

  for (auto const &p : partCfg) {
    auto const r_sq = box_geo.get_mi_vector(pos, p.pos()).norm2();
    if (r_sq < dist_sq) {
      ids.push_back(p.id());
    }
  }

  return ids;
}

std::vector<std::vector<double>>
calc_part_distribution(PartCfg &partCfg, std::vector<int> const &p1_types,
                       std::vector<int> const &p2_types, double r_min,
                       double r_max, int r_bins, bool log_flag, bool int_flag) {

  auto const r_max2 = Utils::sqr(r_max);
  auto const r_min2 = Utils::sqr(r_min);
  auto const start_dist2 = Utils::sqr(r_max + 1.);
  auto const inv_bin_width =
      (log_flag) ? static_cast<double>(r_bins) / std::log(r_max / r_min)
                 : static_cast<double>(r_bins) / (r_max - r_min);

  long cnt = 0;
  double low = 0.0;
  std::vector<double> distribution(r_bins);

  for (auto const &p1 : partCfg) {
    if (Utils::contains(p1_types, p1.type())) {
      auto min_dist2 = start_dist2;
      /* particle loop: p2_types */
      for (auto const &p2 : partCfg) {
        if (p1 != p2) {
          if (Utils::contains(p2_types, p2.type())) {
            auto const act_dist2 =
                box_geo.get_mi_vector(p1.pos(), p2.pos()).norm2();
            if (act_dist2 < min_dist2) {
              min_dist2 = act_dist2;
            }
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
structure_factor(PartCfg &partCfg, std::vector<int> const &p_types, int order) {

  if (order < 1)
    throw std::domain_error("order has to be a strictly positive number");

  auto const order_sq = Utils::sqr(static_cast<std::size_t>(order));
  auto const twoPI_L = 2. * Utils::pi() * box_geo.length_inv()[0];
  std::vector<double> ff(2 * order_sq + 1);

  for (int i = 0; i <= order; i++) {
    for (int j = -order; j <= order; j++) {
      for (int k = -order; k <= order; k++) {
        auto const n = i * i + j * j + k * k;
        if ((static_cast<std::size_t>(n) <= order_sq) && (n >= 1)) {
          double C_sum = 0.0, S_sum = 0.0;
          for (auto const &p : partCfg) {
            if (Utils::contains(p_types, p.type())) {
              auto const qr = twoPI_L * (Utils::Vector3i{{i, j, k}} * p.pos());
              C_sum += cos(qr);
              S_sum += sin(qr);
            }
          }
          ff[2 * n - 2] += C_sum * C_sum + S_sum * S_sum;
          ff[2 * n - 1]++;
        }
      }
    }
  }

  long n_particles = 0l;
  for (auto const &p : partCfg) {
    if (Utils::contains(p_types, p.type())) {
      n_particles++;
    }
  }

  int length = 0;
  for (std::size_t qi = 0; qi < order_sq; qi++) {
    if (ff[2 * qi + 1] != 0) {
      ff[2 * qi] /= static_cast<double>(n_particles) * ff[2 * qi + 1];
      length++;
    }
  }

  std::vector<double> wavevectors(length);
  std::vector<double> intensities(length);

  int cnt = 0;
  for (std::size_t i = 0; i < order_sq; i++) {
    if (ff[2 * i + 1] != 0) {
      wavevectors[cnt] = twoPI_L * sqrt(static_cast<long>(i + 1));
      intensities[cnt] = ff[2 * i];
      cnt++;
    }
  }
  return {std::move(wavevectors), std::move(intensities)};
}
