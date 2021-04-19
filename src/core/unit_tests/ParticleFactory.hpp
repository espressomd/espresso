/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef REACTION_ENSEMBLE_TESTS_PARTICLE_FACTORY_HPP
#define REACTION_ENSEMBLE_TESTS_PARTICLE_FACTORY_HPP

#include "particle_data.hpp"

#include <utils/Vector.hpp>

#include <vector>

/** Fixture to create particles during a test and remove them at the end. */
struct ParticleFactory {
  ParticleFactory() = default;

  ~ParticleFactory() {
    for (auto pid : particle_cache) {
      remove_particle(pid);
    }
  }

  void create_particle(Utils::Vector3d const &pos, int pid = -1,
                       int type = -1) {
    if (pid < 0) {
      pid = get_maximal_particle_id() + 1;
    }
    if (type < 0) {
      type = 0;
    }
    place_particle(pid, pos);
    set_particle_type(pid, type);
    particle_cache.emplace_back(pid);
  }

private:
  std::vector<int> particle_cache;
};

#endif
