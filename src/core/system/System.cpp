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

#include "config/config.hpp"

#include "System.hpp"
#include "System.impl.hpp"

#include "event.hpp"
#include "integrate.hpp"

#include <utils/Vector.hpp>

#include <memory>

namespace System {

static std::shared_ptr<System> instance = std::make_shared<System>();

System::System() {
  box_geo = std::make_shared<BoxGeometry>();
  local_geo = std::make_shared<LocalBox>();
  cell_structure = std::make_shared<CellStructure>(*box_geo);
  nonbonded_ias = std::make_shared<InteractionsNonBonded>();
  comfixed = std::make_shared<ComFixed>();
  force_cap = 0.;
  min_global_cut = INACTIVE_CUTOFF;
}

bool is_system_set() { return instance != nullptr; }

void reset_system() { instance.reset(); }

void set_system(std::shared_ptr<System> new_instance) {
  instance = new_instance;
}

System &get_system() { return *instance; }

Utils::Vector3d System::box() const { return box_geo->length(); }

void System::init() {
#ifdef CUDA
  gpu.init();
#endif
}

void System::set_force_cap(double value) {
  force_cap = value;
  ::recalc_forces = true;
}

void System::set_min_global_cut(double new_value) {
  min_global_cut = new_value;
  on_skin_change();
}

} // namespace System
