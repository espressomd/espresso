/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
    Implementation of \ref statistics_chain.hpp "statistics_chain.hpp".
*/
#include "PartCfg.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "statistics.hpp"
#include <vector>

/** Particles' initial positions (needed for g1(t), g2(t), g3(t)) */
/*@{*/
float *partCoord_g = nullptr, *partCM_g = nullptr;
int n_part_g = 0, n_chains_g = 0;
/*@}*/

/** data for a system consisting of chains. TBRS. */
/*@{*/
int chain_start = 0;
int chain_n_chains = 0;
int chain_length = 0;
/*@}*/

void update_mol_ids_setchains(const ParticleRange &particles) {
  for (auto &p : particles) {
    p.p.mol_id = floor((p.p.identity - chain_start) / (double)chain_length);
  }
}

std::array<double, 4> calc_re(PartCfg &partCfg) {
  double dist = 0.0, dist2 = 0.0, dist4 = 0.0;
  std::array<double, 4> re;

  for (int i = 0; i < chain_n_chains; i++) {
    auto const d =
        partCfg[chain_start + i * chain_length + chain_length - 1].r.p -
        partCfg[chain_start + i * chain_length].r.p;
    auto const norm2 = d.norm2();
    dist += sqrt(norm2);
    dist2 += norm2;
    dist4 += norm2 * norm2;
  }
  auto const tmp = static_cast<double>(chain_n_chains);
  re[0] = dist / tmp;
  re[2] = dist2 / tmp;
  re[1] = sqrt(re[2] - re[0] * re[0]);
  re[3] = sqrt(dist4 / tmp - re[2] * re[2]);
  return re;
}

std::array<double, 4> calc_rg(PartCfg &partCfg) {
  int p;
  double r_G = 0.0, r_G2 = 0.0, r_G4 = 0.0;
  double tmp;
  std::array<double, 4> rg;

  for (int i = 0; i < chain_n_chains; i++) {
    double M = 0.0;
    Utils::Vector3d r_CM{};
    for (int j = 0; j < chain_length; j++) {
      p = chain_start + i * chain_length + j;
      if (partCfg[p].p.is_virtual) {
        throw std::runtime_error(
            "Gyration tensor is not well-defined for chains including virtual "
            "sites. Virtual sites do not have a meaningful mass.");
      }
      r_CM += partCfg[p].r.p * partCfg[p].p.mass;
      M += partCfg[p].p.mass;
    }
    r_CM /= M;
    tmp = 0.0;
    for (int j = 0; j < chain_length; ++j) {
      p = chain_start + i * chain_length + j;
      Utils::Vector3d const d = partCfg[p].r.p - r_CM;
      tmp += d.norm2();
    }
    tmp /= (double)chain_length;
    r_G += sqrt(tmp);
    r_G2 += tmp;
    r_G4 += tmp * tmp;
  }
  tmp = static_cast<double>(chain_n_chains);
  rg[0] = r_G / tmp;
  rg[2] = r_G2 / tmp;
  rg[1] = sqrt(rg[2] - rg[0] * rg[0]);
  rg[3] = sqrt(r_G4 / tmp - rg[2] * rg[2]);
  return rg;
}

std::array<double, 2> calc_rh(PartCfg &partCfg) {
  double r_H = 0.0, r_H2 = 0.0, ri = 0.0, prefac, tmp;
  std::array<double, 2> rh;

  prefac = 0.5 * chain_length *
           (chain_length - 1); /* 1/N^2 is not a normalization factor */
  for (int p = 0; p < chain_n_chains; p++) {
    ri = 0.0;
    for (int i = chain_start + chain_length * p;
         i < chain_start + chain_length * (p + 1); i++) {
      for (int j = i + 1; j < chain_start + chain_length * (p + 1); j++) {
        auto const d = partCfg[i].r.p - partCfg[j].r.p;
        ri += 1.0 / d.norm();
      }
    }
    tmp = prefac / ri;
    r_H += tmp;
    r_H2 += tmp * tmp;
  }
  tmp = static_cast<double>(chain_n_chains);
  rh[0] = r_H / tmp;
  rh[1] = sqrt(r_H2 / tmp - rh[0] * rh[0]);
  return rh;
}
