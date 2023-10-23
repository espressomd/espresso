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
#include "ObservableStat.hpp"

#include "core/BoxGeometry.hpp"
#include "core/analysis/statistics.hpp"
#include "core/analysis/statistics_chain.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/cells.hpp"
#include "core/communication.hpp"
#include "core/dpd.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "script_interface/communication.hpp"

#include <utils/Vector.hpp>
#include <utils/contains.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Analysis {

/** @brief Check if a contiguous range of particle ids exists. */
static void check_topology(CellStructure const &cell_structure, int chain_start,
                           int chain_length, int n_chains) {
  try {
    if (n_chains <= 0) {
      throw std::domain_error("Chain analysis needs at least 1 chain");
    }
    if (chain_length <= 0) {
      throw std::domain_error("Chain analysis needs at least 1 bead per chain");
    }

    auto n_particles_local = 0;
    for (int i = 0; i < n_chains; ++i) {
      for (int j = 0; j < chain_length; ++j) {
        auto const pid = chain_start + i * chain_length + j;
        auto ptr = cell_structure.get_local_particle(pid);
        if (ptr != nullptr and not ptr->is_ghost()) {
          ++n_particles_local;
        }
      }
    }
    auto const n_particles_total = boost::mpi::all_reduce(
        ::comm_cart, n_particles_local, std::plus<int>{});

    if (n_particles_total != n_chains * chain_length) {
      for (int i = 0; i < chain_length * n_chains; ++i) {
        auto const pid = chain_start + i;
        auto ptr = cell_structure.get_local_particle(pid);
        int local_count = 0;
        if (ptr != nullptr and not ptr->is_ghost()) {
          local_count = 1;
        }
        auto const total_count =
            boost::mpi::all_reduce(::comm_cart, local_count, std::plus<int>{});
        if (total_count == 0) {
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
  } catch (...) {
    if (::comm_cart.rank() == 0) {
      throw;
    }
    throw Exception("");
  }
}

void Analysis::check_particle_type(int p_type) const {
  auto const &nonbonded_ias = get_system().nonbonded_ias;
  if (p_type < 0 or p_type > nonbonded_ias->get_max_seen_particle_type()) {
    std::stringstream error_msg;
    error_msg << "Particle type " << p_type << " does not exist";
    throw std::invalid_argument(error_msg.str());
  }
}

Variant Analysis::do_call_method(std::string const &name,
                                 VariantMap const &parameters) {
  if (name == "linear_momentum") {
    auto const local = calc_linear_momentum(
        get_system(), get_value_or<bool>(parameters, "include_particles", true),
        get_value_or<bool>(parameters, "include_lbfluid", true));
    return mpi_reduce_sum(context()->get_comm(), local).as_vector();
  }
  if (name == "particle_energy") {
    auto &system = get_system();
    auto const pid = get_value<int>(parameters, "pid");
    auto const local = system.particle_short_range_energy_contribution(pid);
    return mpi_reduce_sum(context()->get_comm(), local);
  }
#ifdef DIPOLE_FIELD_TRACKING
  if (name == "calc_long_range_fields") {
    get_system().calculate_long_range_fields();
    return {};
  }
#endif
  if (name == "particle_neighbor_pids") {
    get_system().on_observable_calc();
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
#ifdef DPD
  if (name == "dpd_stress") {
    auto const result = dpd_stress(context()->get_comm());
    return result.as_vector();
  }
#endif // DPD
  if (name == "min_dist") {
    auto const p_types1 = get_value<std::vector<int>>(parameters, "p_types1");
    auto const p_types2 = get_value<std::vector<int>>(parameters, "p_types2");
    for (auto const p_type : p_types1) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    for (auto const p_type : p_types2) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    return mindist(get_system(), p_types1, p_types2);
  }
  if (name == "center_of_mass") {
    auto const p_type = get_value<int>(parameters, "p_type");
    context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    auto const local = center_of_mass(get_system(), p_type);
    return mpi_reduce_sum(context()->get_comm(), local).as_vector();
  }
  if (name == "angular_momentum") {
    auto const p_type = get_value<int>(parameters, "p_type");
    context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    auto const local = angular_momentum(get_system(), p_type);
    return mpi_reduce_sum(context()->get_comm(), local).as_vector();
  }
  if (name == "nbhood") {
    auto const pos = get_value<Utils::Vector3d>(parameters, "pos");
    auto const radius = get_value<double>(parameters, "r_catch");
    auto const result = nbhood(get_system(), pos, radius);
    return result;
  }
  if (name == "calc_re") {
    auto const &system = get_system();
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(*system.cell_structure, chain_start, chain_length, n_chains);
    auto const result = calc_re(system, chain_start, chain_length, n_chains);
    return std::vector<double>(result.begin(), result.end());
  }
  if (name == "calc_rg") {
    auto const &system = get_system();
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(*system.cell_structure, chain_start, chain_length, n_chains);
    Variant output;
    context()->parallel_try_catch([&]() {
      auto const result = calc_rg(system, chain_start, chain_length, n_chains);
      output = Variant{std::vector<double>(result.begin(), result.end())};
    });
    return output;
  }
  if (name == "calc_rh") {
    auto const &system = get_system();
    auto const chain_start = get_value<int>(parameters, "chain_start");
    auto const chain_length = get_value<int>(parameters, "chain_length");
    auto const n_chains = get_value<int>(parameters, "number_of_chains");
    check_topology(*system.cell_structure, chain_start, chain_length, n_chains);
    auto const result = calc_rh(system, chain_start, chain_length, n_chains);
    return std::vector<double>(result.begin(), result.end());
  }
  if (name == "gyration_tensor") {
    auto const p_types = get_value<std::vector<int>>(parameters, "p_types");
    for (auto const p_type : p_types) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    auto const mat = gyration_tensor(get_system(), p_types);
    return std::vector<double>(mat.begin(), mat.end());
  }
  if (name == "moment_of_inertia_matrix") {
    auto const p_type = get_value<int>(parameters, "p_type");
    context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    auto const local = moment_of_inertia_matrix(get_system(), p_type);
    return mpi_reduce_sum(context()->get_comm(), local).as_vector();
  }
  if (name == "structure_factor") {
    auto const order = get_value<int>(parameters, "sf_order");
    auto const p_types = get_value<std::vector<int>>(parameters, "sf_types");
    context()->parallel_try_catch([order]() {
      if (order < 1)
        throw std::domain_error("order has to be a strictly positive number");
    });
    for (auto const p_type : p_types) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    auto const result = structure_factor(get_system(), p_types, order);
    return make_vector_of_variants(result);
  }
  if (name == "distribution") {
    auto const &box_l = get_system().box_geo->length();
    auto const r_max_limit =
        0.5 * std::min(std::min(box_l[0], box_l[1]), box_l[2]);
    auto const r_min = get_value_or<double>(parameters, "r_min", 0.);
    auto const r_max = get_value_or<double>(parameters, "r_max", r_max_limit);
    auto const r_bins = get_value_or<int>(parameters, "r_bins", 100);
    auto const log_flag = get_value_or<bool>(parameters, "log_flag", false);
    auto const int_flag = get_value_or<bool>(parameters, "int_flag", false);
    context()->parallel_try_catch([=]() {
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
    });
    auto const p_types1 =
        get_value<std::vector<int>>(parameters, "type_list_a");
    auto const p_types2 =
        get_value<std::vector<int>>(parameters, "type_list_b");
    for (auto const p_type : p_types1) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    for (auto const p_type : p_types2) {
      context()->parallel_try_catch([&]() { check_particle_type(p_type); });
    }
    return make_vector_of_variants(
        calc_part_distribution(get_system(), p_types1, p_types2, r_min, r_max,
                               r_bins, log_flag, int_flag));
  }
  if (name == "calculate_energy") {
    return m_obs_stat->do_call_method("calculate_energy", {});
  }
  if (name == "calculate_scalar_pressure") {
    return m_obs_stat->do_call_method("calculate_scalar_pressure", {});
  }
  if (name == "calculate_pressure_tensor") {
    return m_obs_stat->do_call_method("calculate_pressure_tensor", {});
  }
  return {};
}

} // namespace Analysis
} // namespace ScriptInterface
