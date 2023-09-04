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

#include "electrostatics/solver.hpp"
#include "magnetostatics/solver.hpp"

#include "ek/Solver.hpp"
#include "lb/Solver.hpp"

#include <utils/Vector.hpp>

#include <memory>

struct CellStructure;

namespace System {

class System {
public:
  System();
#ifdef CUDA
  GpuParticleData gpu;
#endif
  ResourceCleanup cleanup_queue;

  Utils::Vector3d box() const;
  void init();

  Coulomb::Solver coulomb;
  Dipoles::Solver dipoles;
  LB::Solver lb;
  EK::Solver ek;
  std::shared_ptr<CellStructure> cell_structure;
};

System &get_system();
void set_system(std::shared_ptr<System> new_instance);
void reset_system();
bool is_system_set();

} // namespace System
