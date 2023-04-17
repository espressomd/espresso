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
 *  Implementation of \ref statistics_chain.hpp "statistics_chain.hpp".
 */

#include "analysis/statistics_chain.hpp"

#include "Particle.hpp"
#include "grid.hpp"
#include "particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

std::array<double, 4> calc_re(int chain_start, int n_chains, int chain_length) {
  double dist = 0.0, dist2 = 0.0, dist4 = 0.0;
  std::array<double, 4> re;

  for (int i = 0; i < n_chains; i++) {
    auto const &p1 =
        get_particle_data(chain_start + i * chain_length + chain_length - 1);
    auto const &p2 = get_particle_data(chain_start + i * chain_length);

    auto const d =
        unfolded_position(p1.pos(), p1.image_box(), box_geo.length()) -
        unfolded_position(p2.pos(), p2.image_box(), box_geo.length());
    auto const norm2 = d.norm2();
    dist += sqrt(norm2);
    dist2 += norm2;
    dist4 += norm2 * norm2;
  }
  auto const tmp = static_cast<double>(n_chains);
  re[0] = dist / tmp;
  re[2] = dist2 / tmp;
  re[1] = (n_chains == 1) ? 0. : std::sqrt(re[2] - Utils::sqr(re[0]));
  re[3] = (n_chains == 1) ? 0. : std::sqrt(dist4 / tmp - Utils::sqr(re[2]));
  return re;
}

std::array<double, 4> calc_rg(int chain_start, int n_chains, int chain_length) {
  double r_G = 0.0, r_G2 = 0.0, r_G4 = 0.0;
  std::array<double, 4> rg;

  for (int i = 0; i < n_chains; i++) {
    double M = 0.0;
    Utils::Vector3d r_CM{};
    for (int j = 0; j < chain_length; j++) {
      auto const &p = get_particle_data(chain_start + i * chain_length + j);

      if (p.is_virtual()) {
        throw std::runtime_error(
            "Gyration tensor is not well-defined for chains including virtual "
            "sites. Virtual sites do not have a meaningful mass.");
      }
      r_CM += unfolded_position(p.pos(), p.image_box(), box_geo.length()) *
              p.mass();
      M += p.mass();
    }
    r_CM /= M;
    double tmp = 0.0;
    for (int j = 0; j < chain_length; ++j) {
      auto const &p = get_particle_data(chain_start + i * chain_length + j);
      Utils::Vector3d const d =
          unfolded_position(p.pos(), p.image_box(), box_geo.length()) - r_CM;
      tmp += d.norm2();
    }
    tmp /= static_cast<double>(chain_length);
    r_G += sqrt(tmp);
    r_G2 += tmp;
    r_G4 += tmp * tmp;
  }
  auto const tmp = static_cast<double>(n_chains);
  rg[0] = r_G / tmp;
  rg[2] = r_G2 / tmp;
  rg[1] = (n_chains == 1) ? 0. : std::sqrt(rg[2] - Utils::sqr(rg[0]));
  rg[3] = (n_chains == 1) ? 0. : std::sqrt(r_G4 / tmp - Utils::sqr(rg[2]));
  return rg;
}

std::array<double, 2> calc_rh(int chain_start, int n_chains, int chain_length) {
  double r_H = 0.0, r_H2 = 0.0;
  std::array<double, 2> rh;

  auto const chain_l = static_cast<double>(chain_length);
  auto const prefac = 0.5 * chain_l * (chain_l - 1.);
  for (int p = 0; p < n_chains; p++) {
    double ri = 0.0;
    for (int i = chain_start + chain_length * p;
         i < chain_start + chain_length * (p + 1); i++) {
      auto const &p1 = get_particle_data(i);
      for (int j = i + 1; j < chain_start + chain_length * (p + 1); j++) {
        auto const &p2 = get_particle_data(j);
        auto const d =
            unfolded_position(p1.pos(), p1.image_box(), box_geo.length()) -
            unfolded_position(p2.pos(), p2.image_box(), box_geo.length());
        ri += 1.0 / d.norm();
      }
    }
    auto const tmp = prefac / ri;
    r_H += tmp;
    r_H2 += tmp * tmp;
  }
  auto const tmp = static_cast<double>(n_chains);
  rh[0] = r_H / tmp;
  rh[1] = (n_chains == 1) ? 0. : std::sqrt(r_H2 / tmp - Utils::sqr(rh[0]));
  return rh;
}
