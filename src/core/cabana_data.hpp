/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include <iostream>

#ifdef CABANA

#include <Cabana_Core.hpp>
#include "custom_verlet_list.hpp"
#include <unordered_map>
#endif

#ifdef CABANA

using data_types = Cabana::MemberTypes<double[3], double[3], int, int, int>;
using memory_space = Kokkos::SharedSpace;
using execution_space = Kokkos::DefaultExecutionSpace;

using ListAlgorithm = Cabana::HalfNeighborTag;
using ListType = Cabana::CustomVerletList<memory_space, ListAlgorithm, Cabana::VerletLayout2D>;

class CabanaData {
    ListType verlet_list;
    std::unordered_map<int, int> id_to_index;

public:
    CabanaData() = default;
    CabanaData(ListType verlet_list, std::unordered_map<int, int> id_to_index)
      : verlet_list(verlet_list), id_to_index(id_to_index) {}

    ListType get_verlet_list() const { return verlet_list; }
    std::unordered_map<int, int> get_id_to_index() const { return id_to_index; }

    ~CabanaData() {
      //std::cout << "Destroying CabanaData" << std::endl;
    };


};
#endif 