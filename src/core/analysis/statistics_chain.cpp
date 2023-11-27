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

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

/**
 * @brief Gather particle properties (or any derived quantities) on MPI rank 0.
 * The data is collected into an unordered map.
 */
template <typename T> struct GatherParticleTraits {
private:
  std::vector<int> buffer_pid{};
  std::vector<T> buffer_obs{};

public:
  virtual T kernel(Particle const &) const = 0;

  void fetch(CellStructure const &cell_structure, int pid) {
    auto const ptr = cell_structure.get_local_particle(pid);
    if (ptr != nullptr and not ptr->is_ghost()) {
      buffer_pid.emplace_back(pid);
      buffer_obs.emplace_back(kernel(*ptr));
    }
  }

  auto join() {
    std::unordered_map<int, T> map{};
    Utils::Mpi::gather_buffer(buffer_pid, ::comm_cart, 0);
    Utils::Mpi::gather_buffer(buffer_obs, ::comm_cart, 0);
    if (::comm_cart.rank() == 0) {
      for (std::size_t i = 0u; i < buffer_pid.size(); ++i) {
        map[buffer_pid[i]] = buffer_obs[i];
      }
    }
    buffer_pid.clear();
    buffer_obs.clear();
    return map;
  }
};

struct GatherPos : public GatherParticleTraits<Utils::Vector3d> {
  BoxGeometry const &m_box_geo;
  GatherPos(BoxGeometry const &box_geo) : m_box_geo{box_geo} {}
  Utils::Vector3d kernel(Particle const &p) const override {
    return m_box_geo.unfolded_position(p.pos(), p.image_box());
  }
};

struct GatherCom : public GatherParticleTraits<Utils::Vector3d> {
  BoxGeometry const &m_box_geo;
  GatherCom(BoxGeometry const &box_geo) : m_box_geo{box_geo} {}
  Utils::Vector3d kernel(Particle const &p) const override {
    return m_box_geo.unfolded_position(p.pos(), p.image_box()) * p.mass();
  }
};

struct GatherMass : public GatherParticleTraits<double> {
  double kernel(Particle const &p) const override { return p.mass(); }
};

std::array<double, 4> calc_re(System::System const &system, int chain_start,
                              int chain_length, int n_chains) {
  auto const &cell_structure = *system.cell_structure;
  GatherPos prefetch{*system.box_geo};
  double dist = 0.0, dist2 = 0.0, dist4 = 0.0;
  std::array<double, 4> re{};

  for (int i = 0; i < n_chains; i++) {
    auto const pid2 = chain_start + i * chain_length;
    auto const pid1 = pid2 + chain_length - 1;
    prefetch.fetch(cell_structure, pid1);
    prefetch.fetch(cell_structure, pid2);
  }
  auto const map = prefetch.join();
  if (::comm_cart.rank() == 0) {
    for (int i = 0; i < n_chains; i++) {
      auto const pid2 = chain_start + i * chain_length;
      auto const pid1 = pid2 + chain_length - 1;
      auto const norm2 = (map.at(pid1) - map.at(pid2)).norm2();
      dist += sqrt(norm2);
      dist2 += norm2;
      dist4 += norm2 * norm2;
    }
    auto const tmp = static_cast<double>(n_chains);
    re[0] = dist / tmp;
    re[2] = dist2 / tmp;
    re[1] = (n_chains == 1) ? 0. : std::sqrt(re[2] - Utils::sqr(re[0]));
    re[3] = (n_chains == 1) ? 0. : std::sqrt(dist4 / tmp - Utils::sqr(re[2]));
  }
  return re;
}

std::array<double, 4> calc_rg(System::System const &system, int chain_start,
                              int chain_length, int n_chains) {
  auto const &cell_structure = *system.cell_structure;
  GatherPos prefetch_pos{*system.box_geo};
  GatherCom prefetch_com{*system.box_geo};
  GatherMass prefetch_mass{};
  double r_G = 0.0, r_G2 = 0.0, r_G4 = 0.0;
  std::array<double, 4> rg{};

  auto has_virtual = false;
  for (int i = 0; i < n_chains * chain_length; ++i) {
    auto const pid = chain_start + i;
    auto const ptr = cell_structure.get_local_particle(pid);
    if (ptr != nullptr and not ptr->is_ghost() and ptr->is_virtual()) {
      has_virtual = true;
      break;
    }
  }
  if (boost::mpi::all_reduce(::comm_cart, has_virtual, std::logical_or<>{})) {
    throw std::runtime_error(
        "Center of mass is not well-defined for chains including virtual "
        "sites. Virtual sites do not have a meaningful mass.");
  }

  for (int i = 0; i < n_chains * chain_length; ++i) {
    auto const pid = chain_start + i;
    prefetch_com.fetch(cell_structure, pid);
    prefetch_pos.fetch(cell_structure, pid);
    prefetch_mass.fetch(cell_structure, pid);
  }

  auto const map_pos = prefetch_pos.join();
  auto const map_com = prefetch_com.join();
  auto const map_mass = prefetch_mass.join();
  if (::comm_cart.rank() == 0) {
    for (int i = 0; i < n_chains; i++) {
      double M = 0.0;
      Utils::Vector3d r_CM{};
      for (int j = 0; j < chain_length; j++) {
        auto const pid = chain_start + i * chain_length + j;
        r_CM += map_com.at(pid);
        M += map_mass.at(pid);
      }
      r_CM /= M;
      double tmp = 0.0;
      for (int j = 0; j < chain_length; ++j) {
        auto const pid = chain_start + i * chain_length + j;
        auto const d = map_pos.at(pid) - r_CM;
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
  }
  return rg;
}

std::array<double, 2> calc_rh(System::System const &system, int chain_start,
                              int chain_length, int n_chains) {
  auto const &cell_structure = *system.cell_structure;
  GatherPos prefetch{*system.box_geo};
  double r_H = 0.0, r_H2 = 0.0;
  std::array<double, 2> rh{};

  auto const chain_l = static_cast<double>(chain_length);
  auto const prefac = 0.5 * chain_l * (chain_l - 1.);
  for (int p = 0; p < n_chains; p++) {
    for (int i = chain_start + chain_length * p;
         i < chain_start + chain_length * (p + 1); i++) {
      prefetch.fetch(cell_structure, i);
      for (int j = i + 1; j < chain_start + chain_length * (p + 1); j++) {
        prefetch.fetch(cell_structure, j);
      }
    }
  }
  auto const map = prefetch.join();
  if (::comm_cart.rank() == 0) {
    for (int p = 0; p < n_chains; p++) {
      double ri = 0.0;
      for (int i = chain_start + chain_length * p;
           i < chain_start + chain_length * (p + 1); i++) {
        for (int j = i + 1; j < chain_start + chain_length * (p + 1); j++) {
          ri += 1.0 / (map.at(i) - map.at(j)).norm();
        }
      }
      auto const tmp = prefac / ri;
      r_H += tmp;
      r_H2 += tmp * tmp;
    }
    auto const tmp = static_cast<double>(n_chains);
    rh[0] = r_H / tmp;
    rh[1] = (n_chains == 1) ? 0. : std::sqrt(r_H2 / tmp - Utils::sqr(rh[0]));
  }
  return rh;
}
