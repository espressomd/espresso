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

#include "System.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <memory>

namespace System {

static std::shared_ptr<System> instance;

void reset_system() { instance.reset(); }

void set_system(std::shared_ptr<System> new_instance) {
  instance = new_instance;
}

System &get_system() { return *instance; }

Utils::Vector3d System::box() const { return ::box_geo.length(); }

} // namespace System
