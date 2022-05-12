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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_CELL_SYSTEM_CELL_SYSTEM_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_CELL_SYSTEM_CELL_SYSTEM_HPP

#include "script_interface/ScriptInterface.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/cells.hpp"
#include "core/communication.hpp"
#include "core/event.hpp"
#include "core/grid.hpp"
#include "core/integrate.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/particle_node.hpp"
#include "core/tuning.hpp"

#include <utils/Vector.hpp>
#include <utils/as_const.hpp>

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

namespace {
Variant pack_vector(Utils::Vector3i const &vec) {
  return std::vector<int>(vec.begin(), vec.end());
}
auto const &get_regular_decomposition() {
  return dynamic_cast<RegularDecomposition const &>(
      Utils::as_const(cell_structure).decomposition());
}

auto const &get_hybrid_decomposition() {
  return dynamic_cast<HybridDecomposition const &>(
      Utils::as_const(cell_structure).decomposition());
}
} // namespace

class CellSystem : public AutoParameters<CellSystem> {
  std::unordered_map<CellStructureType, std::string> const cs_type_to_name = {
      {CellStructureType::CELL_STRUCTURE_REGULAR, "regular_decomposition"},
      {CellStructureType::CELL_STRUCTURE_NSQUARE, "n_square"},
      {CellStructureType::CELL_STRUCTURE_HYBRID, "hybrid_decomposition"},
  };

  std::unordered_map<std::string, CellStructureType> const cs_name_to_type = {
      {"regular_decomposition", CellStructureType::CELL_STRUCTURE_REGULAR},
      {"n_square", CellStructureType::CELL_STRUCTURE_NSQUARE},
      {"hybrid_decomposition", CellStructureType::CELL_STRUCTURE_HYBRID},
  };

public:
  CellSystem() {
    add_parameters({
        {"use_verlet_lists", cell_structure.use_verlet_list},
        {"node_grid",
         [this](Variant const &v) {
           context()->parallel_try_catch([&v]() {
             auto const error_msg = std::string("Parameter 'node_grid'");
             auto const vec = get_value<std::vector<int>>(v);
             if (vec.size() != 3ul) {
               throw std::invalid_argument(error_msg + " must be 3 ints");
             }
             auto const new_node_grid = Utils::Vector3i{vec.begin(), vec.end()};
             auto const n_nodes_old = Utils::product(::node_grid);
             auto const n_nodes_new = Utils::product(new_node_grid);
             if (n_nodes_new != n_nodes_old) {
               std::stringstream reason;
               reason << ": MPI world size " << n_nodes_old << " incompatible "
                      << "with new node grid [" << new_node_grid << "]";
               throw std::invalid_argument(error_msg + reason.str());
             }
             ::node_grid = new_node_grid;
             on_nodegrid_change();
           });
         },
         []() { return pack_vector(::node_grid); }},
        {"skin",
         [this](Variant const &v) {
           auto const new_skin = get_value<double>(v);
           if (new_skin < 0.) {
             if (context()->is_head_node()) {
               throw std::domain_error("Parameter 'skin' must be >= 0");
             }
             throw Exception("");
           }
           mpi_set_skin_local(new_skin);
         },
         []() { return ::skin; }},
        {"decomposition_type", AutoParameter::read_only,
         [this]() {
           return cs_type_to_name.at(cell_structure.decomposition_type());
         }},
        {"n_square_types", AutoParameter::read_only,
         []() {
           if (cell_structure.decomposition_type() !=
               CellStructureType::CELL_STRUCTURE_HYBRID) {
             return Variant{none};
           }
           auto const hd = get_hybrid_decomposition();
           auto const ns_types = hd.get_n_square_types();
           return Variant{std::vector<int>(ns_types.begin(), ns_types.end())};
         }},
        {"cutoff_regular", AutoParameter::read_only,
         []() {
           if (cell_structure.decomposition_type() !=
               CellStructureType::CELL_STRUCTURE_HYBRID) {
             return Variant{none};
           }
           auto const hd = get_hybrid_decomposition();
           return Variant{hd.get_cutoff_regular()};
         }},
        {"max_cut_nonbonded", AutoParameter::read_only,
         maximal_cutoff_nonbonded},
        {"max_cut_bonded", AutoParameter::read_only, maximal_cutoff_bonded},
        {"interaction_range", AutoParameter::read_only, interaction_range},
    });
  }

  void do_construct(VariantMap const &params) override {
    // the core cell structure has a default constructor; params is empty
    // during default construction for the python class, and not empty
    // during checkpointing
    if (params.count("decomposition_type")) {
      auto const cs_name = get_value<std::string>(params, "decomposition_type");
      auto const cs_type = cs_name_to_type.at(cs_name);
      initialize(cs_type, params);
      do_set_parameter("skin", params.at("skin"));
      do_set_parameter("node_grid", params.at("node_grid"));
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
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
      auto const cs_type = cell_structure.decomposition_type();
      if (cs_type == CellStructureType::CELL_STRUCTURE_REGULAR) {
        auto const rd = get_regular_decomposition();
        state["cell_grid"] = pack_vector(rd.cell_grid);
        state["cell_size"] = Variant{rd.cell_size};
      } else if (cs_type == CellStructureType::CELL_STRUCTURE_HYBRID) {
        auto const hd = get_hybrid_decomposition();
        state["cell_grid"] = pack_vector(hd.get_cell_grid());
        state["cell_size"] = Variant{hd.get_cell_size()};
        mpi_resort_particles(true); // needed to get correct particle counts
        auto const n_particles = hybrid_parts_per_decomposition(hd);
        state["parts_per_decomposition"] =
            Variant{std::unordered_map<std::string, Variant>{
                {"regular", n_particles.first},
                {"n_square", n_particles.second}}};
      }
      state["verlet_reuse"] = get_verlet_reuse();
      state["n_nodes"] = ::n_nodes;
      return state;
    }
    if (name == "get_pairs") {
      std::vector<Variant> out;
      context()->parallel_try_catch([&params, &out]() {
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
        std::transform(pair_list.begin(), pair_list.end(),
                       std::back_inserter(out),
                       [](std::pair<int, int> const &pair) {
                         return std::vector<int>{pair.first, pair.second};
                       });
      });
      return out;
    }
    if (name == "get_neighbors") {
      std::vector<std::vector<int>> neighbors_global;
      context()->parallel_try_catch([&neighbors_global, &params]() {
        auto const dist = get_value<double>(params, "distance");
        auto const pid = get_value<int>(params, "pid");
        auto const ret = mpi_get_short_range_neighbors_local(pid, dist, true);
        std::vector<int> neighbors_local;
        if (ret) {
          neighbors_local = *ret;
        }
        boost::mpi::gather(comm_cart, neighbors_local, neighbors_global, 0);
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
      std::vector<Variant> out;
      auto const pair_list = non_bonded_loop_trace();
      std::transform(pair_list.begin(), pair_list.end(),
                     std::back_inserter(out), [](PairInfo const &pair) {
                       return std::vector<Variant>{pair.id1,   pair.id2,
                                                   pair.pos1,  pair.pos2,
                                                   pair.vec21, pair.node};
                     });
      return out;
    }
    if (name == "tune_skin") {
      if (context()->is_head_node()) {
        tune_skin(get_value<double>(params, "min_skin"),
                  get_value<double>(params, "max_skin"),
                  get_value<double>(params, "tol"),
                  get_value<int>(params, "int_steps"),
                  get_value_or<bool>(params, "adjust_max_skin", false));
      }
      return ::skin;
    }
    return {};
  }

private:
  /**
   * @brief Resort the particles.
   *
   * This function resorts the particles on the nodes.
   *
   * @param global_flag   If true a global resort is done, if false particles
   *                      are only exchanges between neighbors.
   * @return The number of particles on the nodes after the resort.
   */
  std::vector<int> mpi_resort_particles(bool global_flag) const {
    cell_structure.resort_particles(global_flag, box_geo);
    if (context()->is_head_node()) {
      clear_particle_node();
    }
    auto const size = static_cast<int>(cell_structure.local_particles().size());
    std::vector<int> n_part_per_node;
    boost::mpi::gather(comm_cart, size, n_part_per_node, 0);
    return n_part_per_node;
  }

  void initialize(CellStructureType const &cs_type,
                  VariantMap const &params) const {
    auto const verlet = get_value_or<bool>(params, "use_verlet_lists", true);
    cell_structure.use_verlet_list = verlet;
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
};

} // namespace CellSystem
} // namespace ScriptInterface

#endif
