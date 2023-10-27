/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include <system/System.hpp>

#include <utils/Vector.hpp>

#include <memory>

namespace System {
class System;
}

/** Manager for a stand-alone ESPResSo system.
 *  The system is default-initialized, MPI-ready and has no script interface.
 */
class EspressoSystemStandAlone {
public:
  EspressoSystemStandAlone(int argc, char **argv);
  void set_box_l(Utils::Vector3d const &box_l) const;
  void set_node_grid(Utils::Vector3i const &node_grid) const;
  void set_time_step(double time_step) const;
  void set_skin(double new_skin) const;
  auto get_handle() { return m_instance; }

private:
  bool head_node;
  std::shared_ptr<System::System> m_instance;
};
