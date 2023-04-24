/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "CellSystem.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/cell_system/HybridDecomposition.hpp"
#include "core/cell_system/RegularDecomposition.hpp"
#include "core/cells.hpp"
#include "core/event.hpp"
#include "core/grid.hpp"
#include "core/integrate.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/particle_node.hpp"
#include "core/tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace CellSystem {

static auto const &get_regular_decomposition() {
  return dynamic_cast<RegularDecomposition const &>(
      std::as_const(::cell_structure).decomposition());
}

static auto const &get_hybrid_decomposition() {
  return dynamic_cast<HybridDecomposition const &>(
      std::as_const(::cell_structure).decomposition());
}

CellSystem::CellSystem() {
  add_parameters({
      {"use_verlet_lists", ::cell_structure.use_verlet_list},
      {"node_grid",
       [this](Variant const &v) {
         context()->parallel_try_catch([&v]() {
           auto const error_msg = std::string("Parameter 'node_grid'");
           auto const old_node_grid = ::node_grid;
           auto const new_node_grid = get_value<Utils::Vector3i>(v);
           auto const n_nodes_old = Utils::product(old_node_grid);
           auto const n_nodes_new = Utils::product(new_node_grid);
           if (n_nodes_new != n_nodes_old) {
             std::stringstream reason;
             reason << ": MPI world size " << n_nodes_old << " incompatible "
                    << "with new node grid [" << new_node_grid << "]";
             throw std::invalid_argument(error_msg + reason.str());
           }
           try {
             set_node_grid(new_node_grid);
           } catch (...) {
             set_node_grid(old_node_grid);
             throw;
           }
         });
       },
       []() { return ::node_grid; }},
      {"skin",
       [this](Variant const &v) {
         auto const new_skin = get_value<double>(v);
         if (new_skin < 0.) {
           if (context()->is_head_node()) {
             throw std::domain_error("Parameter 'skin' must be >= 0");
           }
           throw Exception("");
         }
         ::set_skin(new_skin);
       },
       []() { return ::skin; }},
      {"decomposition_type", AutoParameter::read_only,
       [this]() {
         return cs_type_to_name.at(::cell_structure.decomposition_type());
       }},
      {"n_square_types", AutoParameter::read_only,
       []() {
         if (::cell_structure.decomposition_type() !=
             CellStructureType::CELL_STRUCTURE_HYBRID) {
           return Variant{none};
         }
         auto const hd = get_hybrid_decomposition();
         auto const ns_types = hd.get_n_square_types();
         return Variant{std::vector<int>(ns_types.begin(), ns_types.end())};
       }},
      {"cutoff_regular", AutoParameter::read_only,
       []() {
         if (::cell_structure.decomposition_type() !=
             CellStructureType::CELL_STRUCTURE_HYBRID) {
           return Variant{none};
         }
         auto const hd = get_hybrid_decomposition();
         return Variant{hd.get_cutoff_regular()};
       }},
      {"max_cut_nonbonded", AutoParameter::read_only, maximal_cutoff_nonbonded},
      {"max_cut_bonded", AutoParameter::read_only, maximal_cutoff_bonded},
      {"interaction_range", AutoParameter::read_only, interaction_range},
  });
}

Variant CellSystem::do_call_method(std::string const &name,
                                   VariantMap const &params) {
  if (name == "initialize") {
    auto const cs_name = get_value<std::string>(params, "name");
    auto const cs_type = cs_name_to_type.at(cs_name);
    initialize(cs_type, params);
    return {};
  }
  if (name == "resort") {
    auto const global_flag = get_value_or<bool>(params, "global_flag", true);
    return mpi_resort_particles(global_flag);
  }
  if (name == "get_state") {
    auto state = get_parameters();
    auto const cs_type = ::cell_structure.decomposition_type();
    if (cs_type == CellStructureType::CELL_STRUCTURE_REGULAR) {
      auto const rd = get_regular_decomposition();
      state["cell_grid"] = Variant{rd.cell_grid};
      state["cell_size"] = Variant{rd.cell_size};
    } else if (cs_type == CellStructureType::CELL_STRUCTURE_HYBRID) {
      auto const hd = get_hybrid_decomposition();
      state["cell_grid"] = Variant{hd.get_cell_grid()};
      state["cell_size"] = Variant{hd.get_cell_size()};
      mpi_resort_particles(true); // needed to get correct particle counts
      state["parts_per_decomposition"] =
          Variant{std::unordered_map<std::string, Variant>{
              {"regular", hd.count_particles_in_regular()},
              {"n_square", hd.count_particles_in_n_square()}}};
    }
    state["verlet_reuse"] = get_verlet_reuse();
    state["n_nodes"] = context()->get_comm().size();
    return state;
  }
  if (name == "get_pairs") {
    on_observable_calc();
    std::vector<Variant> out;
    context()->parallel_try_catch([this, &params, &out]() {
      std::vector<std::pair<int, int>> pair_list;
      auto const distance = get_value<double>(params, "distance");
      if (boost::get<std::string>(&params.at("types")) != nullptr) {
        auto const key = get_value<std::string>(params, "types");
        if (key != "all") {
          throw std::invalid_argument("Unknown argument types='" + key + "'");
        }
        pair_list = get_pairs(distance);
      } else {
        auto const types = get_value<std::vector<int>>(params, "types");
        pair_list = get_pairs_of_types(distance, types);
      }
      Utils::Mpi::gather_buffer(pair_list, context()->get_comm());
      std::transform(pair_list.begin(), pair_list.end(),
                     std::back_inserter(out),
                     [](std::pair<int, int> const &pair) {
                       return std::vector<int>{pair.first, pair.second};
                     });
    });
    return out;
  }
  if (name == "get_neighbors") {
    on_observable_calc();
    std::vector<std::vector<int>> neighbors_global;
    context()->parallel_try_catch([this, &neighbors_global, &params]() {
      auto const dist = get_value<double>(params, "distance");
      auto const pid = get_value<int>(params, "pid");
      auto const ret = get_short_range_neighbors(pid, dist);
      std::vector<int> neighbors_local;
      if (ret) {
        neighbors_local = *ret;
      }
      boost::mpi::gather(context()->get_comm(), neighbors_local,
                         neighbors_global, 0);
    });
    std::vector<int> neighbors;
    for (auto const &neighbors_local : neighbors_global) {
      if (not neighbors_local.empty()) {
        neighbors = neighbors_local;
        break;
      }
    }
    return neighbors;
  }
  if (name == "non_bonded_loop_trace") {
    on_observable_calc();
    std::vector<Variant> out;
    auto pair_list = non_bonded_loop_trace(context()->get_comm().rank());
    Utils::Mpi::gather_buffer(pair_list, context()->get_comm());
    std::transform(pair_list.begin(), pair_list.end(), std::back_inserter(out),
                   [](PairInfo const &pair) {
                     return std::vector<Variant>{pair.id1,   pair.id2,
                                                 pair.pos1,  pair.pos2,
                                                 pair.vec21, pair.node};
                   });
    return out;
  }
  if (name == "tune_skin") {
    tune_skin(get_value<double>(params, "min_skin"),
              get_value<double>(params, "max_skin"),
              get_value<double>(params, "tol"),
              get_value<int>(params, "int_steps"),
              get_value_or<bool>(params, "adjust_max_skin", false));
    return ::skin;
  }
  if (name == "get_max_range") {
    return ::cell_structure.max_range();
  }
  return {};
}

std::vector<int> CellSystem::mpi_resort_particles(bool global_flag) const {
  ::cell_structure.resort_particles(global_flag, box_geo);
  clear_particle_node();
  auto const size = static_cast<int>(::cell_structure.local_particles().size());
  std::vector<int> n_part_per_node;
  boost::mpi::gather(context()->get_comm(), size, n_part_per_node, 0);
  return n_part_per_node;
}

void CellSystem::initialize(CellStructureType const &cs_type,
                            VariantMap const &params) const {
  auto const verlet = get_value_or<bool>(params, "use_verlet_lists", true);
  ::cell_structure.use_verlet_list = verlet;
  if (cs_type == CellStructureType::CELL_STRUCTURE_HYBRID) {
    auto const cutoff_regular = get_value<double>(params, "cutoff_regular");
    auto const ns_types =
        get_value_or<std::vector<int>>(params, "n_square_types", {});
    auto n_square_types = std::set<int>{ns_types.begin(), ns_types.end()};
    set_hybrid_decomposition(std::move(n_square_types), cutoff_regular);
  } else {
    cells_re_init(cs_type);
  }
}

} // namespace CellSystem
} // namespace ScriptInterface
