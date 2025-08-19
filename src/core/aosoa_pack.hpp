/*
 * Copyright (C) 2025 The ESPResSo project
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

#include <config/config.hpp>

#ifdef SHARED_MEMORY_PARALLELISM

#include "cell_system/CellStructure.hpp"

#include <Cabana_Core.hpp>

struct CellStructure::AoSoA_pack {
  CellStructure::AoSoAType::member_slice_type<0> position;
  CellStructure::AoSoAType::member_slice_type<1> charge;
  CellStructure::AoSoAType::member_slice_type<2> id;
  CellStructure::AoSoAType::member_slice_type<3> type;

  AoSoA_pack() = default;

  AoSoA_pack(CellStructure::AoSoAType &aosoa)
      : position(Cabana::slice<0>(aosoa)), charge(Cabana::slice<1>(aosoa)),
        id(Cabana::slice<2>(aosoa)), type(Cabana::slice<3>(aosoa)) {}
};

#endif // SHARED_MEMORY_PARALLELISM
