/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "h5md_specification.hpp"
#include "h5md_core.hpp"
#include "hdf5.h"

#include <utility>

namespace Writer {
namespace H5md {

H5MD_Specification::H5MD_Specification(unsigned int fields) {
  auto const add_time_series = [this](Dataset &&dataset, bool link = true) {
    auto const group = dataset.group;
    m_datasets.push_back(std::move(dataset));
    m_datasets.push_back({group, "step", 1, H5T_NATIVE_INT, 1, link});
    m_datasets.push_back({group, "time", 1, H5T_NATIVE_DOUBLE, 1, link});
  };

  if (fields & H5MD_OUT_BOX_L) {
    add_time_series(
        {"particles/atoms/box/edges", "value", 2, H5T_NATIVE_DOUBLE, 3, false});
  }
  if (fields & H5MD_OUT_LE_OFF) {
    add_time_series({"particles/atoms/lees_edwards/offset", "value", 2,
                     H5T_NATIVE_DOUBLE, 1, false});
  }
  if (fields & H5MD_OUT_LE_DIR) {
    add_time_series({"particles/atoms/lees_edwards/direction", "value", 2,
                     H5T_NATIVE_INT, 1, false});
  }
  if (fields & H5MD_OUT_LE_NORMAL) {
    add_time_series({"particles/atoms/lees_edwards/normal", "value", 2,
                     H5T_NATIVE_INT, 1, false});
  }
  if (fields & H5MD_OUT_MASS) {
    add_time_series(
        {"particles/atoms/mass", "value", 2, H5T_NATIVE_DOUBLE, 1, false});
  }
  if (fields & H5MD_OUT_CHARGE) {
    add_time_series(
        {"particles/atoms/charge", "value", 2, H5T_NATIVE_DOUBLE, 1, false});
  }
  add_time_series({"particles/atoms/id", "value", 2, H5T_NATIVE_INT, 1, false},
                  false);
  if (fields & H5MD_OUT_TYPE) {
    add_time_series(
        {"particles/atoms/species", "value", 2, H5T_NATIVE_INT, 1, false});
  }
  if (fields & H5MD_OUT_POS) {
    add_time_series(
        {"particles/atoms/position", "value", 3, H5T_NATIVE_DOUBLE, 3, false});
  }
  if (fields & H5MD_OUT_VEL) {
    add_time_series(
        {"particles/atoms/velocity", "value", 3, H5T_NATIVE_DOUBLE, 3, false});
  }
  if (fields & H5MD_OUT_FORCE) {
    add_time_series(
        {"particles/atoms/force", "value", 3, H5T_NATIVE_DOUBLE, 3, false});
  }
  if (fields & H5MD_OUT_IMG) {
    add_time_series(
        {"particles/atoms/image", "value", 3, H5T_NATIVE_INT, 3, false});
  }
  if (fields & H5MD_OUT_BONDS) {
    add_time_series(
        {"connectivity/atoms", "value", 3, H5T_NATIVE_INT, 2, false});
  }
}

} // namespace H5md
} // namespace Writer
