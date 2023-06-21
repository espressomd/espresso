/*
 * Copyright (C) 2013-2022 The ESPResSo project
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

#include "Analysis.hpp"

#include "core/analysis/statistics.hpp"
#include "core/analysis/statistics_chain.hpp"
#include "core/cells.hpp"
#include "core/dpd.hpp"
#include "core/energy.hpp"
#include "core/event.hpp"
#include "core/grid.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/partCfg_global.hpp"
#include "core/particle_node.hpp"

#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Analysis {

/** @brief Check if a contiguous range of particle ids exists. */
static void check_topology(int chain_start, int chain_length, int n_chains) {
  if (n_chains <= 0) {
    throw std::domain_error("Chain analysis needs at least 1 chain");
  }
  if (chain_length <= 0) {
    throw std::domain_error("Chain analysis needs at least 1 bead per chain");
  }
  for (int i = 0; i < chain_length * n_chains; ++i) {
    auto const pid = chain_start + i;
    if (not particle_exists(pid)) {
      std::stringstream error_msg;
      error_msg << "Particle with id " << pid << " does not exist; "
                << "cannot perform analysis on the range chain_start="
                << chain_start << ", number_of_chains=" << n_chains
                << ", chain_length=" << chain_length << ". "
                << "Please provide a contiguous range of particle ids.";
      throw std::runtime_error(error_msg.str());
    }
  }
}

/** @brief Check if a particle type exists. */
static void check_particle_type(int p_type) {
  if (p_type < 0 or p_type >= ::max_seen_particle_type) {
    std::stringstream error_msg;
    error_msg << "Particle type " << p_type << " does not exist";
    throw std::invalid_argument(error_msg.str());
  }
}

Variant Analysis::do_call_method(std::string const &name,
                                 VariantMap const &parameters) {
  if (name == "linear_momentum") {
    auto const local = calc_linear_momentum(
        get_value_or<bool>(parameters, "include_particles", true),
        get_value_or<bool>(parameters, "include_lbfluid", true));
    return mpi_reduce_sum(context()->get_comm(), local).as_vector();
  }
  if (name == "particle_energy") {
    auto const pid = get_value<int>(parameters, "pid");
    auto const local = particle_short_range_energy_contribution(pid);
    return mpi_reduce_sum(context()->get_comm(), local);
  }
#ifdef DIPOLE_FIELD_TRACKING
  if (name == "calc_long_range_fields") {
    calc_long_range_fields();
    return {};
  }
#endif
  if (name == "particle_neighbor_pids") {
    on_observable_calc();
    std::unordered_map<int, std::vector<int>> dict;
    context()->parallel_try_catch([&]() {
      auto neighbor_pids = get_neighbor_pids();
      Utils::Mpi::gather_buffer(neighbor_pids, context()->get_comm());
      std::for_each(neighbor_pids.begin(), neighbor_pids.end(),
                    [&dict](NeighborPIDs const &neighbor_pid) {
                      dict[neighbor_pid.pid] = neighbor_pid.neighbor_pids;
                    });
    });
    return make_unordered_map_of_variants(dict);
  }
  if (not context()->is_head_node()) {
    return {};
  }
  if (name == "min_dist") {
    auto const p_types1 = get_value<std::vector<int>>(parameters, "p_types1");
    auto const p_types2 = get_value<std::vector<int>>(parameters, "p_types2");
    for (auto const p_type : p_types1) {
      check_particle_type(p_type);
    }
    for (auto const p_type : p_types2) {
      check_particle_type(p_type);
    }
    return mindist(partCfg(), p_types1, p_types2);
  }
  if (name == "center_of_mass") {
    auto const p_type = get_value<int>(parameters, "p_type");
    check_particle_type(p_type);
    return center_of_mass(partCfg(), p_type).as_vector();
  }
  if (name == "angular_momentum") {
    auto const p_type = get_value<int>(parameters, "p_type");
    auto const result = angular_momentum(partCfg(), p_type);
    return result.as_vector();
  }
  if (name == "nbhood") {
    auto const pos = get_value<Utils::Vector3d>(parameters, "pos");
    auto const radius = get_value<double>(parameters, "r_catch");
    auto const result = nbhood(partCfg(), pos, radius);
    return result;
  }
#ifdef DPD
  if (name == "dpd_stress") {
    auto const result = dpd_stress();
    return result.as_vector();
  }
#endif // DPD
  if (name == "calc_re") {
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(chain_start, chain_length, n_chains);
    auto const result = calc_re(chain_start, n_chains, chain_length);
    return std::vector<double>(result.begin(), result.end());
  }
  if (name == "calc_rg") {
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(chain_start, chain_length, n_chains);
    auto const result = calc_rg(chain_start, n_chains, chain_length);
    return std::vector<double>(result.begin(), result.end());
  }
  if (name == "calc_rh") {
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(chain_start, chain_length, n_chains);
    auto const result = calc_rh(chain_start, n_chains, chain_length);
    return std::vector<double>(result.begin(), result.end());
  }
  if (name == "gyration_tensor") {
    auto const p_types = get_value<std::vector<int>>(parameters, "p_types");
    for (auto const p_type : p_types) {
      check_particle_type(p_type);
    }
    std::vector<Utils::Vector3d> positions{};
    for (const auto &p : partCfg()) {
      if (Utils::contains(p_types, p.type())) {
        positions.push_back(p.pos());
      }
    }
    auto const com =
        std::accumulate(positions.begin(), positions.end(), Utils::Vector3d{}) /
        static_cast<double>(positions.size());
    // compute covariance matrix
    Utils::Vector9d mat{};
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i > j) {
          mat[i * 3 + j] = mat[j * 3 + i];
        } else {
          mat[i * 3 + j] = std::accumulate(
              positions.begin(), positions.end(), 0.,
              [i, j, &com](double acc, Utils::Vector3d const &pos) {
                return acc + (pos[i] - com[i]) * (pos[j] - com[j]);
              });
        }
      }
    }
    mat /= static_cast<double>(positions.size());
    return std::vector<double>(mat.begin(), mat.end());
  }
  if (name == "moment_of_inertia_matrix") {
    auto const p_type = get_value<int>(parameters, "p_type");
    check_particle_type(p_type);
    auto const result = moment_of_inertia_matrix(partCfg(), p_type);
    return result.as_vector();
  }
  if (name == "structure_factor") {
    auto const order = get_value<int>(parameters, "sf_order");
    auto const p_types = get_value<std::vector<int>>(parameters, "sf_types");
    for (auto const p_type : p_types) {
      check_particle_type(p_type);
    }
    auto const result = structure_factor(partCfg(), p_types, order);
    return make_vector_of_variants(result);
  }
  if (name == "distribution") {
    auto const r_max_limit =
        0.5 * std::min(std::min(::box_geo.length()[0], ::box_geo.length()[1]),
                       ::box_geo.length()[2]);
    auto const r_min = get_value_or<double>(parameters, "r_min", 0.);
    auto const r_max = get_value_or<double>(parameters, "r_max", r_max_limit);
    auto const r_bins = get_value_or<int>(parameters, "r_bins", 100);
    auto const log_flag = get_value_or<bool>(parameters, "log_flag", false);
    auto const int_flag = get_value_or<bool>(parameters, "int_flag", false);
    if (log_flag and r_min <= 0.) {
      throw std::domain_error("Parameter 'r_min' must be > 0");
    }
    if (r_min < 0.) {
      throw std::domain_error("Parameter 'r_min' must be >= 0");
    }
    if (r_min >= r_max) {
      throw std::domain_error("Parameter 'r_max' must be > 'r_min'");
    }
    if (r_max > r_max_limit) {
      throw std::domain_error("Parameter 'r_max' must be <= box_l / 2");
    }
    if (r_bins <= 0) {
      throw std::domain_error("Parameter 'r_bins' must be >= 1");
    }
    auto const p_types1 =
        get_value<std::vector<int>>(parameters, "type_list_a");
    auto const p_types2 =
        get_value<std::vector<int>>(parameters, "type_list_b");
    for (auto const p_type : p_types1) {
      check_particle_type(p_type);
    }
    for (auto const p_type : p_types2) {
      check_particle_type(p_type);
    }
    return make_vector_of_variants(
        calc_part_distribution(partCfg(), p_types1, p_types2, r_min, r_max,
                               r_bins, log_flag, int_flag));
  }
  return {};
}

} // namespace Analysis
} // namespace ScriptInterface
