/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#include "ReactionAlgorithm.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/particle_data/ParticleHandle.hpp"
#include "script_interface/system/System.hpp"

#include "core/Observable_stat.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/system/System.hpp"

#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/serialization.hpp>

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace ReactionMethods {

Variant ReactionAlgorithm::do_call_method(std::string const &name,
                                          VariantMap const &params) {
  if (name == "count_number_of_particles_per_type") {
    auto const &cs = m_cell_system->get_cell_structure();
    auto const types = get_value<std::vector<int>>(params, "types");
    std::vector<int> local_numbers;
    std::vector<int> global_numbers(types.size());
    context()->parallel_try_catch([&]() {
      for (auto const &type : types) {
        if (type < 0) {
          throw std::runtime_error("Types may not be negative");
        }
        int counter = 0;
        for (auto const &p : cs.local_particles()) {
          if (p.type() == type) {
            counter++;
          }
        }
        local_numbers.emplace_back(counter);
      }
    });
    boost::mpi::reduce(context()->get_comm(), local_numbers, global_numbers,
                       std::plus<>(), 0);
    return global_numbers;
  }
  if (name == "single_update") {
    auto const pid = get_value<int>(params, "pid");
    auto const properties = get_value<VariantMap>(params, "properties");
    m_particle_modifier->set_pid(pid);
    auto const &cs = m_cell_system->get_cell_structure();
    auto p = cs.get_local_particle(pid);
    if (p.has_value() and p->is_ghost()) {
      p.reset();
    }
    int old_type = -1;
    if (context()->is_head_node()) {
      if (p) {
        old_type = p->type();
      } else {
        context()->get_comm().recv(boost::mpi::any_source, 42, old_type);
      }
    } else if (p) {
      context()->get_comm().send(0, 42, p->type());
    }
    for (auto const &[param_name, value] :
         std::map<std::string, Variant>(properties.begin(), properties.end())) {
      m_particle_modifier->do_set_parameter(param_name, value);
    }
    return old_type;
  }
  if (name == "batch_update") {
    auto const pids = get_value<std::vector<int>>(params, "pids");
    auto const properties = get_value<VariantMap>(params, "properties");
    for (int pid : pids) {
      m_particle_modifier->set_pid(pid);
      for (auto const &[param_name, value] : std::map<std::string, Variant>(
               properties.begin(), properties.end())) {
        m_particle_modifier->do_set_parameter(param_name, value);
      }
    }
    return {};
  }
  if (name == "delete_particle") {
    m_particle_modifier->set_pid(get_value<int>(params, "pid"));
    m_particle_modifier->ParticleHandle::do_call_method("remove_particle", {});
    return {};
  }
  if (name == "delete_particles") {
    auto const pids = get_value<std::vector<int>>(params, "pids");
    for (int pid : pids) {
      m_particle_modifier->set_pid(pid);
      m_particle_modifier->ParticleHandle::do_call_method("remove_particle",
                                                          {});
    }
    return {};
  }
  if (context()->is_head_node()) {
    throw std::runtime_error("unknown method '" + name + "'");
  }
  return {};
}

void ReactionAlgorithm::do_construct(VariantMap const &params) {
  m_system = get_value<std::shared_ptr<System::System>>(params, "system");
  m_cell_system = get_value<std::shared_ptr<CellSystem::CellSystem>>(
      m_system->get_parameter("cell_system"));
  m_particle_modifier = get_value<std::shared_ptr<Particles::ParticleModifier>>(
      params, "particle_modifier");
}

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */
