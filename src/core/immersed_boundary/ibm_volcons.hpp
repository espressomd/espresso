/*
 * Copyright (C) 2010-2021 The ESPResSo project
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

#ifndef IBM_VOLCONS_H
#define IBM_VOLCONS_H

#include <utils/Vector.hpp>

/** Parameters for IBM volume conservation bond */
struct IBMVolCons {
  /** ID of the large soft particle to which this node belongs */
  int softID;
  /** Reference volume */
  double volRef;
  /** Spring constant for volume force */
  double kappaV;

  double cutoff() const { return 0.; }

  static constexpr int num = 0;

  IBMVolCons(int softID, double kappaV);

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &softID;
    ar &volRef;
    ar &kappaV;
  }
};

#endif
