/*
 * Copyright (C) 2014-2022 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#include "GpuParticleData.hpp"
#include "ResourceCleanup.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <memory>

namespace System {

class System {
public:
#ifdef CUDA
  GpuParticleData gpu;
#endif
  ResourceCleanup cleanup_queue;

  Utils::Vector3d box() const;
  void init() {
#ifdef CUDA
    gpu.init();
#endif
  }
};

System &get_system();
void set_system(std::shared_ptr<System> new_instance);
void reset_system();

} // namespace System
