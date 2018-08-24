/*
  Copyright (C) 2017 The ESPResSo project

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
#ifndef CORE_COMFIXED_HPP
#define CORE_COMFIXED_HPP

#include "Vector.hpp"
#include "utils/keys.hpp"

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>

#include <unordered_map>
#include <vector>

template <typename ParticleRange> class ComFixed {
public:
  using TypeIndex = std::unordered_map<int, int>;

private:
  TypeIndex m_type_index;

  std::vector<Vector3d> local_type_forces(ParticleRange &particles) const {
    std::vector<Vector3d> ret(m_type_index.size(), Vector3d{});

    for (auto const &p : particles) {
      /* Check if type is of interest */
      auto it = m_type_index.find(p.p.type);
      if (it != m_type_index.end()) {
        ret[it->second] += Vector3d{p.f.f};
      }
    }

    return ret;
  }

  std::vector<double> local_type_masses(ParticleRange &particles) const {
    std::vector<double> ret(m_type_index.size(), 0.);

    for (auto const &p : particles) {
      /* Check if type is of interest */
      auto it = m_type_index.find(p.p.type);
      if (it != m_type_index.end()) {
        ret[it->second] += p.p.mass;
      }
    }

    return ret;
  }

public:
  template <typename Container> void set_fixed_types(Container const &c) {
    m_type_index.clear();

    int i = 0;
    for (auto const &t : c) {
      m_type_index[t] = i++;
    }
  }
  std::vector<int> get_fixed_types() const { return Utils::keys(m_type_index); }

  void apply(boost::mpi::communicator const &comm,
             ParticleRange &particles) const {
    /* Bail out early if there is nothing to do. */
    if (m_type_index.empty())
      return;

    auto const local_forces = local_type_forces(particles);
    auto const local_masses = local_type_masses(particles);

    /* Total forces and masses of the types. */
    std::vector<Vector3d> forces(m_type_index.size(), Vector3d{});
    std::vector<double> masses(m_type_index.size(), 0.0);

    /* Add contributions from all nodes and redistribute them to all. */
    boost::mpi::all_reduce(comm, local_forces.data(), local_forces.size(),
                           forces.data(), std::plus<Vector3d>{});
    boost::mpi::all_reduce(comm, local_masses.data(), local_masses.size(),
                           masses.data(), std::plus<double>{});

    for (auto &p : particles) {
      /* Check if type is of interest */
      auto it = m_type_index.find(p.p.type);
      if (it != m_type_index.end()) {
        auto const mass_frac = p.p.mass / masses[it->second];
        auto const &type_force = forces[it->second];
        for (int i = 0; i < 3; i++) {
          p.f.f[i] -= mass_frac * type_force[i];
        }
      }
    }
  }
};

#endif
