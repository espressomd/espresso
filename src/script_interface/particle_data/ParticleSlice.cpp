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

#include "ParticleSlice.hpp"
#include "ParticleHandle.hpp"

#include "core/system/System.hpp"

#include "core/bonds.hpp"
#include "core/particle_node.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "script_interface/Context.hpp"
#include "script_interface/Exception.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"
#include "utils/Vector.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

void set_particles_bonds(
    std::vector<int> const pids,
    std::vector<std::vector<int>> const all_bonds_ids,
    std::vector<std::vector<std::vector<int>>> const all_bonds_partner_ids,
    ::CellStructure &cell_structure, Context *context,
    std::shared_ptr<::System::System> system) {
  for (std::size_t i = 0; i < pids.size(); i += 1) {
    auto const pid = pids[i];
    auto const bonds_ids = all_bonds_ids[i];
    auto const bonds_partner_ids = all_bonds_partner_ids[i];

    // Remove old bonds
    auto p_ptr = get_real_particle(context->get_comm(), pid, cell_structure);
    if (p_ptr != nullptr) {
      p_ptr->bonds().clear();
    }
    // Add new bonds
    for (std::size_t j = 0; j < bonds_ids.size(); j += 1) {
      std::vector<int> particle_ids = {pid};
      std::ranges::copy(bonds_partner_ids[j], std::back_inserter(particle_ids));
      ::add_bond(*system, bonds_ids[j], particle_ids);
      system->on_particle_change();
    }
  }
}

#ifdef ESPRESSO_EXCLUSIONS
void set_particles_exclusions(
    std::vector<int> const pids,
    std::vector<std::vector<int>> const exclusion_lists,
    ::CellStructure &cell_structure, Context *context,
    std::shared_ptr<::System::System> system) {
  for (std::size_t i = 0; i < pids.size(); i += 1) {
    auto const pid = pids[i];
    auto const exclusion_list = exclusion_lists[i];
    context->parallel_try_catch([&]() {
      for (auto const excluded_pid : exclusion_list) {
        particle_exclusion_sanity_checks(pid, excluded_pid, cell_structure,
                                         context);
      }
    });
    auto p_ptr = get_real_particle(context->get_comm(), pid, cell_structure);
    if (p_ptr != nullptr) {
      auto const &p = *p_ptr;
      // Remove all excluded ids of this particle
      for (auto const old_excluded_pid : p.exclusions()) {
        local_remove_exclusion(pid, old_excluded_pid, cell_structure);
      }
      // Add new excluded ids for this particle
      for (auto const excluded_pid : exclusion_list) {
        if (!p.has_exclusion(excluded_pid)) {
          local_add_exclusion(pid, excluded_pid, cell_structure);
        }
      }
    }
  } // Particle id loop
  system->on_particle_change();
}
#endif // ESPRESSO_EXCLUSIONS

void set_particles_positions(std::vector<int> const pids,
                             std::vector<Utils::Vector3d> const positions) {
  for (std::size_t i = 0; i < pids.size(); i += 1) {
    auto const pid = pids[i];
    auto const pos = positions[i];
    particle_checks(pid, pos);
    set_particle_pos(pid, pos);
  } // Particle id loop
}

void set_particles_types(std::vector<int> const pids,
                         std::vector<int> const types,
                         boost::mpi::communicator const &comm,
                         CellStructure &cell_structure,
                         std::shared_ptr<::System::System> system) {
  for (std::size_t i = 0; i < pids.size(); i += 1) {
    auto const pid = pids[i];
    auto p_ptr = get_real_particle(comm, pid, cell_structure);
    if (p_ptr != nullptr) {
      auto const old_type = p_ptr->type();
      auto const &new_type = types[i];
      if (new_type < 0) {
        throw std::domain_error(error_msg("type", "must be an integer >= 0"));
      }
      system->nonbonded_ias->make_particle_type_exist(new_type);
      on_particle_type_change(pid, old_type, new_type);
      p_ptr->type() = new_type;
    }
  } // Particle id loop
}

#ifdef ESPRESSO_ELECTROSTATICS
void set_particles_charges(
    std::vector<int> const pids,
    std::variant<std::vector<int>, std::vector<double>> const &charges,
    boost::mpi::communicator const &comm, CellStructure &cell_structure,
    std::shared_ptr<::System::System> system) {
  std::visit(
      [&](auto const &concrete_charges) {
        for (std::size_t i = 0; i < pids.size(); i += 1) {
          auto const pid = pids[i];
          auto p_ptr = get_real_particle(comm, pid, cell_structure);
          if (p_ptr != nullptr) {
            p_ptr->q() = concrete_charges[i];
          }
        } // Particle id loop
        system->on_particle_charge_change();
      },
      charges);
}
#endif // ESPRESSO_ELECTROSTATICS

void ParticleSlice::do_construct(VariantMap const &params) {
  if (params.contains("__cell_structure")) {
    auto so = get_value<std::shared_ptr<CellSystem::CellSystem>>(
        params, "__cell_structure");
    so->configure(*this);
    m_cell_structure = so;
  }
  if (params.contains("__bonded_ias")) {
    m_bonded_ias = get_value<std::shared_ptr<Interactions::BondedInteractions>>(
        params, "__bonded_ias");
  }
  m_id_selection = get_value<std::vector<int>>(params, "id_selection");
  m_chunk_size = get_value_or<int>(params, "prefetch_chunk_size", 10000);
  if (not context()->is_head_node()) {
    return;
  }
  for (auto const pid : m_id_selection) {
    if (not particle_exists(pid)) {
      throw std::out_of_range("Particle does not exist: " +
                              std::to_string(pid));
    }
  }
}

Variant ParticleSlice::do_call_method(std::string const &name,
                                      VariantMap const &params) {
  if (name == "set_param_parallel") {
    auto const param_name = get_value<std::string>(params, "name");
    if (not params.contains("values")) {
      if (param_name == "bonds") {
        if (not params.contains("all_bonds_ids")) {
          throw Exception("Parameter 'all_bonds_ids' is missing.");
        }
        if (not params.contains("all_bonds_partner_ids")) {
          throw Exception("Parameter 'all_bonds_partner_ids' is missing.");
        }
      } else {
        throw Exception("Parameter '" + param_name + "' is missing.");
      }
    }
    // Handle parameters with special setters
    if (m_special_parameters.contains(param_name)) {
      context()->parallel_try_catch([&]() {
        if (param_name == "pos") {
          set_particles_positions(
              m_id_selection,
              get_value<std::vector<Utils::Vector3d>>(params, "values"));
        }

        else if (param_name == "type") {
          set_particles_types(
              m_id_selection, get_value<std::vector<int>>(params, "values"),
              context()->get_comm(), *get_cell_structure(), get_system());
        }
#ifdef ESPRESSO_ELECTROSTATICS
        else if (param_name == "q") {
          std::variant<std::vector<int>, std::vector<double>> charges;
          if (std::holds_alternative<std::vector<int>>(params.at("values"))) {
            charges = get_value<std::vector<int>>(params, "values");
          } else {
            charges = get_value<std::vector<double>>(params, "values");
          }
          set_particles_charges(m_id_selection, charges, context()->get_comm(),
                                *get_cell_structure(), get_system());

        }
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_EXCLUSIONS
        else if (param_name == "exclusions") {
          auto const excluded_pids =
              get_value<std::vector<std::vector<int>>>(params, "values");
          set_particles_exclusions(
              m_id_selection,
              get_value<std::vector<std::vector<int>>>(params, "values"),
              *get_cell_structure(), context(), get_system());
        }
#endif // ESPRESSO_EXCLUSIONS
        else if (param_name == "bonds") {
          set_particles_bonds(
              m_id_selection,
              get_value<std::vector<std::vector<int>>>(params, "all_bonds_ids"),
              get_value<std::vector<std::vector<std::vector<int>>>>(
                  params, "all_bonds_partner_ids"),
              *get_cell_structure(), context(), get_system());
        }
      });
    }
    // Handle generic parameters
    else {
      context()->parallel_try_catch([&]() {
        std::visit(
            [&](auto &&vals) {
              SetParticleParametersVisitor{}(m_id_selection, param_name, vals,
                                             context(), m_cell_structure.lock(),
                                             m_bonded_ias.lock());
            },
            params.at("values"));
      });
    }
    return {};
  }
  if (not context()->is_head_node()) {
    return {};
  }
  if (name == "prefetch_particle_data") {
    auto p_ids = get_value<std::vector<int>>(params, "chunk");
    prefetch_particle_data(p_ids);
    return {};
  }
  if (name == "get_particle") {
    return context()->make_shared(
        "Particles::ParticleHandle",
        {{"id", get_value<int>(params, "p_id")},
         {"__cell_structure", m_cell_structure.lock()},
         {"__bonded_ias", m_bonded_ias.lock()}});
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
