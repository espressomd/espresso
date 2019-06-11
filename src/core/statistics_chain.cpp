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
    Implementation of \ref statistics_chain.hpp "statistics_chain.hpp".
*/
#include "PartCfg.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "statistics.hpp"

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

void update_mol_ids_setchains() {
  for (auto &p : local_cells.particles()) {
    p.p.mol_id = floor((p.p.identity - chain_start) / (double)chain_length);
  }
}

void calc_re(PartCfg &partCfg, double **_re) {
  int i;
  double dx, dy, dz;
  double dist = 0.0, dist2 = 0.0, dist4 = 0.0;
  double *re = nullptr, tmp;
  *_re = re = Utils::realloc(re, 4 * sizeof(double));

  for (i = 0; i < chain_n_chains; i++) {
    dx = partCfg[chain_start + i * chain_length + chain_length - 1].r.p[0] -
         partCfg[chain_start + i * chain_length].r.p[0];
    dy = partCfg[chain_start + i * chain_length + chain_length - 1].r.p[1] -
         partCfg[chain_start + i * chain_length].r.p[1];
    dz = partCfg[chain_start + i * chain_length + chain_length - 1].r.p[2] -
         partCfg[chain_start + i * chain_length].r.p[2];
    tmp = (Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz));
    dist += sqrt(tmp);
    dist2 += tmp;
    dist4 += tmp * tmp;
  }
  tmp = (double)chain_n_chains;
  re[0] = dist / tmp;
  re[2] = dist2 / tmp;
  re[1] = sqrt(re[2] - re[0] * re[0]);
  re[3] = sqrt(dist4 / tmp - re[2] * re[2]);
}

void calc_rg(PartCfg &partCfg, double **_rg) {
  int i, j, p;
  double dx, dy, dz, r_CM_x, r_CM_y, r_CM_z;
  double r_G = 0.0, r_G2 = 0.0, r_G4 = 0.0;
  double *rg = nullptr, IdoubMPC, tmp;
  double M;
  *_rg = rg = Utils::realloc(rg, 4 * sizeof(double));

  for (i = 0; i < chain_n_chains; i++) {
    M = 0.0;
    r_CM_x = r_CM_y = r_CM_z = 0.0;
    IdoubMPC = 1. / (double)chain_length;
    for (j = 0; j < chain_length; j++) {
      p = chain_start + i * chain_length + j;
      r_CM_x += partCfg[p].r.p[0] * (partCfg[p]).p.mass;
      r_CM_y += partCfg[p].r.p[1] * (partCfg[p]).p.mass;
      r_CM_z += partCfg[p].r.p[2] * (partCfg[p]).p.mass;
      M += (partCfg[p]).p.mass;
    }
    r_CM_x /= M;
    r_CM_y /= M;
    r_CM_z /= M;
    tmp = 0.0;
    for (j = 0; j < chain_length; ++j) {
      p = chain_start + i * chain_length + j;
      dx = partCfg[p].r.p[0] - r_CM_x;
      dy = partCfg[p].r.p[1] - r_CM_y;
      dz = partCfg[p].r.p[2] - r_CM_z;
      tmp += (Utils::sqr(dx) + Utils::sqr(dy) + Utils::sqr(dz));
    }
    tmp *= IdoubMPC;
    r_G += sqrt(tmp);
    r_G2 += tmp;
    r_G4 += tmp * tmp;
  }
  tmp = (double)chain_n_chains;
  rg[0] = r_G / tmp;
  rg[2] = r_G2 / tmp;
  rg[1] = sqrt(rg[2] - rg[0] * rg[0]);
  rg[3] = sqrt(r_G4 / tmp - rg[2] * rg[2]);
}

void calc_rh(PartCfg &partCfg, double **_rh) {
  int i, j, p;
  double dx, dy, dz, r_H = 0.0, r_H2 = 0.0, *rh = nullptr, ri = 0.0, prefac,
                     tmp;
  *_rh = rh = Utils::realloc(rh, 2 * sizeof(double));

  prefac = 0.5 * chain_length *
           (chain_length - 1); /* 1/N^2 is not a normalization factor */
  for (p = 0; p < chain_n_chains; p++) {
    ri = 0.0;
    for (i = chain_start + chain_length * p;
         i < chain_start + chain_length * (p + 1); i++) {
      for (j = i + 1; j < chain_start + chain_length * (p + 1); j++) {
        dx = partCfg[i].r.p[0] - partCfg[j].r.p[0];
        dx *= dx;
        dy = partCfg[i].r.p[1] - partCfg[j].r.p[1];
        dy *= dy;
        dz = partCfg[i].r.p[2] - partCfg[j].r.p[2];
        dz *= dz;
        ri += 1.0 / sqrt(dx + dy + dz);
      }
    }
    tmp = prefac / ri;
    r_H += tmp;
    r_H2 += tmp * tmp;
  }
  tmp = (double)chain_n_chains;
  rh[0] = r_H / tmp;
  rh[1] = sqrt(r_H2 / tmp - rh[0] * rh[0]);
}
