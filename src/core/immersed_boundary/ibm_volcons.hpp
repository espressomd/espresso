/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <boost/serialization/access.hpp>

#include <cassert>
#include <stdexcept>
#include <vector>

/** Parameters for IBM volume conservation bond */
struct IBMVolCons {
  /** ID of the large soft particle to which this node belongs */
  unsigned int softID;
  /** Reference volume */
  double volRef;
  /** Spring constant for volume force */
  double kappaV;

  double cutoff() const { return 0.; }

  static constexpr int num = 0;

  IBMVolCons(int softID, double kappaV) {
    if (softID < 0) {
      throw std::domain_error("IBMVolCons parameter 'softID' has to be >= 0");
    }
    this->softID = static_cast<unsigned int>(softID);
    this->kappaV = kappaV;
    // NOTE: We cannot compute the reference volume here because not all
    // interactions are setup and thus we do not know which triangles belong to
    // this softID. Calculate it later in the init function of
    // \ref ImmersedBoundaries::init_volume_conservation()
    volRef = 0.;
    m_volumes = nullptr;
  }

  double get_current_volume() const {
    double volume = 0.;
    if (m_volumes) {
      assert(static_cast<std::size_t>(softID) < m_volumes->size());
      volume = (*m_volumes)[softID];
    }
    return volume;
  }

  void set_volumes_view(std::vector<double> const &volumes) {
    m_volumes = &volumes;
  }
  void unset_volumes_view() { m_volumes = nullptr; }

private:
  std::vector<double> const *m_volumes;
};
