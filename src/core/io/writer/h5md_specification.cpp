/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#include "hdf5.h"

#include <array>

namespace Writer {
namespace H5md {

std::array<H5MD_Specification::Dataset, 39> H5MD_Specification::DATASETS = {{
    {"particles/atoms/box/edges", "value", 2, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/box/edges", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/box/edges", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/lees_edwards/offset", "value", 2, H5T_NATIVE_DOUBLE, 1,
     false},
    {"particles/atoms/lees_edwards/offset", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/lees_edwards/offset", "time", 1, H5T_NATIVE_DOUBLE, 1,
     true},
    {"particles/atoms/lees_edwards/direction", "value", 2, H5T_NATIVE_INT, 1,
     false},
    {"particles/atoms/lees_edwards/direction", "step", 1, H5T_NATIVE_INT, 1,
     true},
    {"particles/atoms/lees_edwards/direction", "time", 1, H5T_NATIVE_DOUBLE, 1,
     true},
    {"particles/atoms/lees_edwards/normal", "value", 2, H5T_NATIVE_INT, 1,
     false},
    {"particles/atoms/lees_edwards/normal", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/lees_edwards/normal", "time", 1, H5T_NATIVE_DOUBLE, 1,
     true},
    {"particles/atoms/mass", "value", 2, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/mass", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/mass", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/charge", "value", 2, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/charge", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/charge", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/id", "value", 2, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/id", "step", 1, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/id", "time", 1, H5T_NATIVE_DOUBLE, 1, false},
    {"particles/atoms/species", "value", 2, H5T_NATIVE_INT, 1, false},
    {"particles/atoms/species", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/species", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/position", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/position", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/position", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/velocity", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/velocity", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/velocity", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/force", "value", 3, H5T_NATIVE_DOUBLE, 3, false},
    {"particles/atoms/force", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/force", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"particles/atoms/image", "value", 3, H5T_NATIVE_INT, 3, false},
    {"particles/atoms/image", "step", 1, H5T_NATIVE_INT, 1, true},
    {"particles/atoms/image", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
    {"connectivity/atoms", "value", 3, H5T_NATIVE_INT, 2, false},
    {"connectivity/atoms", "step", 1, H5T_NATIVE_INT, 1, true},
    {"connectivity/atoms", "time", 1, H5T_NATIVE_DOUBLE, 1, true},
}};
}
} // namespace Writer
