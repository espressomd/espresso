/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "core/bonds.hpp"
#include "core/nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "core/particle_node.hpp"
#include "core/system/System.hpp"

#include "script_interface/Context.hpp"
#include "script_interface/Exception.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>

#include <algorithm>
#include <functional>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace ScriptInterface {
namespace Particles {

static void set_particles_bonds(
    std::vector<int> const &pids,
    std::vector<std::vector<int>> const &all_bonds_ids,
    std::vector<std::vector<std::vector<int>>> const &all_bonds_partner_ids,
    ::CellStructure &cell_structure, ::System::System &system) {
  for (std::size_t i = 0; i < pids.size(); ++i) {
    auto const pid = pids[i];
    auto const bonds_ids = all_bonds_ids[i];
    auto const bonds_partner_ids = all_bonds_partner_ids[i];
    // Remove old bonds
    auto p = cell_structure.get_local_particle(pid);
    if (p != nullptr and not p->is_ghost()) {
      p->bonds().clear();
    }
    // Add new bonds
    for (std::size_t j = 0; j < bonds_ids.size(); ++j) {
      std::vector<int> particle_ids = {pid};
      std::ranges::copy(bonds_partner_ids[j], std::back_inserter(particle_ids));
      ::add_bond(system, bonds_ids[j], particle_ids);
      system.on_particle_change();
    }
  }
}

#ifdef ESPRESSO_EXCLUSIONS
static void
set_particles_exclusions(std::vector<int> const &pids,
                         std::vector<std::vector<int>> const &exclusion_lists,
                         boost::mpi::communicator const &comm,
                         ::CellStructure &cell_structure,
                         ::System::System &system) {
  for (std::size_t i = 0; i < pids.size(); ++i) {
    auto const pid = pids[i];
    auto const &exclusion_list = exclusion_lists[i];
    for (auto const excluded_pid : exclusion_list) { // collective communication
      particle_exclusion_sanity_checks(pid, excluded_pid, cell_structure, comm);
    }
    auto p = cell_structure.get_local_particle(pid);
    if (p != nullptr and not p->is_ghost()) {
      // Remove all excluded ids of this particle
      for (auto const old_excluded_pid : p->exclusions()) {
        local_remove_exclusion(pid, old_excluded_pid, cell_structure);
      }
      // Add new excluded ids for this particle
      for (auto const excluded_pid : exclusion_list) {
        if (not p->has_exclusion(excluded_pid)) {
          local_add_exclusion(pid, excluded_pid, cell_structure);
        }
      }
    }
  }
  system.on_particle_change();
}
#endif // ESPRESSO_EXCLUSIONS

static void
set_particles_positions(std::vector<int> const &pids,
                        std::vector<Utils::Vector3d> const &positions) {
  for (std::size_t i = 0; i < pids.size(); ++i) {
    auto const pid = pids[i];
    auto const pos = positions[i];
    particle_checks(pid, pos);
    set_particle_pos(pid, pos);
  }
}

static void set_particles_types(std::vector<int> const &pids,
                                std::vector<int> const &types,
                                CellStructure &cell_structure,
                                ::System::System &system) {
  for (std::size_t i = 0; i < pids.size(); ++i) {
    auto const pid = pids[i];
    auto p = cell_structure.get_local_particle(pid);
    if (p != nullptr and not p->is_ghost()) {
      auto const old_type = p->type();
      auto const &new_type = types[i];
      if (new_type < 0) {
        throw std::domain_error(error_msg("type", "must be an integer >= 0"));
      }
      system.nonbonded_ias->make_particle_type_exist(new_type);
      on_particle_type_change(pid, old_type, new_type);
      p->type() = new_type;
    }
  }
}

#ifdef ESPRESSO_ELECTROSTATICS
static void set_particles_charges(std::vector<int> const &pids,
                                  std::vector<double> const &charges,
                                  CellStructure &cell_structure,
                                  ::System::System &system) {
  for (std::size_t i = 0; i < pids.size(); ++i) {
    auto const pid = pids[i];
    auto p = cell_structure.get_local_particle(pid);
    if (p != nullptr and not p->is_ghost()) {
      p->q() = charges[i];
    }
  }
  system.on_particle_charge_change();
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
      context()->parallel_try_catch([&]() {
        if (param_name == "bonds") {
          if (not params.contains("all_bonds_ids")) {
            throw Exception("Parameter 'all_bonds_ids' is missing");
          }
          if (not params.contains("all_bonds_partner_ids")) {
            throw Exception("Parameter 'all_bonds_partner_ids' is missing");
          }
        } else {
          throw Exception("Parameter 'values' is missing");
        }
      });
    }
    // Handle parameters with special setters
    if (m_special_parameters.contains(param_name)) {
      context()->parallel_try_catch([&]() {
        if (param_name == "pos") {
          set_particles_positions(
              m_id_selection,
              get_value<std::vector<Utils::Vector3d>>(params, "values"));
        } else if (param_name == "type") {
          set_particles_types(m_id_selection,
                              get_value<std::vector<int>>(params, "values"),
                              *get_cell_structure(), *get_system());
        }
#ifdef ESPRESSO_ELECTROSTATICS
        else if (param_name == "q") {
          std::vector<double> charges;
          if (is_type<std::vector<int>>(params.at("values"))) {
            auto tmp = get_value<std::vector<int>>(params, "values");
            charges = std::vector<double>(tmp.begin(), tmp.end());
          } else {
            charges = get_value<std::vector<double>>(params, "values");
          }
          set_particles_charges(m_id_selection, charges, *get_cell_structure(),
                                *get_system());

        }
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_EXCLUSIONS
        else if (param_name == "exclusions") {
          auto const excluded_pids =
              get_value<std::vector<std::vector<int>>>(params, "values");
          set_particles_exclusions(
              m_id_selection,
              get_value<std::vector<std::vector<int>>>(params, "values"),
              context()->get_comm(), *get_cell_structure(), *get_system());
        }
#endif // ESPRESSO_EXCLUSIONS
        else if (param_name == "bonds") {
          set_particles_bonds(
              m_id_selection,
              get_value<std::vector<std::vector<int>>>(params, "all_bonds_ids"),
              get_value<std::vector<std::vector<std::vector<int>>>>(
                  params, "all_bonds_partner_ids"),
              *get_cell_structure(), *get_system());
        }
      });
    } else {
      // Handle generic parameters
      context()->parallel_try_catch([&]() {
        std::visit(
            [&](auto &&vals) {
              set_from_vector_like(m_id_selection, param_name, vals, context(),
                                   m_cell_structure.lock(),
                                   m_bonded_ias.lock());
            },
            params.at("values"));
      });
    }
    return {};
  }
  if (name == "get_param_parallel") {
    auto const param_name = get_value<std::string>(params, "name");

    // handle special optimized properties
    if (param_name == "type") {
      auto const getter{[](Particle const &p) { return p.type(); }};
      return get_particles_properties<int>(m_id_selection, getter, context(),
                                           *get_cell_structure());
    }
    if (param_name == "q") {
      auto const getter{[](Particle const &p) { return p.q(); }};
      return get_particles_properties<double>(m_id_selection, getter, context(),
                                              *get_cell_structure());
    }
    if (param_name == "pos") {
      auto const &box_geo = *get_system()->box_geo;
      auto const getter = [&box_geo](Particle const &p) {
        return box_geo.unfolded_position(p.pos(), p.image_box());
      };
      return get_particles_properties<Utils::Vector3d>(
          m_id_selection, getter, context(), *get_cell_structure());
    }
    if (param_name == "pos_folded") {
      auto const &box_geo = *get_system()->box_geo;
      auto const getter = [&box_geo](Particle const &p) {
        return box_geo.folded_position(p.pos());
      };
      return get_particles_properties<Utils::Vector3d>(
          m_id_selection, getter, context(), *get_cell_structure());
    }

    // handle all other particle properties using expensive MPI reductions
    if (!context()->is_head_node()) {
      return {};
    }
    std::vector<Variant> result;
    result.reserve(m_id_selection.size());
    auto so = std::dynamic_pointer_cast<ParticleModifier>(
        context()->make_shared("Particles::ParticleModifier",
                               {{"id", -1},
                                {"__cell_structure", m_cell_structure.lock()},
                                {"__bonded_ias", m_bonded_ias.lock()}}));
    for (int pid : m_id_selection) {
      so->set_pid(pid);
      result.emplace_back(so->get_parameter(param_name));
    }
    return result;
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
