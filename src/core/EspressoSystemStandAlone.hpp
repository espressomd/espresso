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
#ifndef ESPRESSO_SYSTEM_STAND_ALONE_HPP
#define ESPRESSO_SYSTEM_STAND_ALONE_HPP

#include <utils/Vector.hpp>

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

private:
  bool head_node;
};

#endif
