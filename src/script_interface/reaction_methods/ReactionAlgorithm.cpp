/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "SingleReaction.hpp"

#include "config/config.hpp"

#include "core/cells.hpp"
#include "core/communication.hpp"
#include "core/energy.hpp"
#include "core/event.hpp"
#include "core/particle_node.hpp"
#include "core/reaction_methods/ReactionAlgorithm.hpp"
#include "core/reaction_methods/utils.hpp"

#include <utils/contains.hpp>

#include <boost/serialization/serialization.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace ScriptInterface {
namespace ReactionMethods {

ReactionAlgorithm::ReactionAlgorithm() {
  add_parameters(
      {{"reactions", AutoParameter::read_only,
        [this]() {
          std::vector<Variant> out;
          for (auto const &e : m_reactions) {
            out.emplace_back(e);
          }
          return out;
        }},
       {"kT", AutoParameter::read_only, [this]() { return RE()->get_kT(); }},
       {"search_algorithm",
        [this](Variant const &v) {
          m_exclusion->do_set_parameter("search_algorithm", v);
        },
        [this]() { return m_exclusion->get_parameter("search_algorithm"); }},
       {"particle_inside_exclusion_range_touched",
        [this](Variant const &v) {
          RE()->particle_inside_exclusion_range_touched = get_value<bool>(v);
        },
        [this]() { return RE()->particle_inside_exclusion_range_touched; }},
       {"default_charges", AutoParameter::read_only,
        [this]() {
          return make_unordered_map_of_variants(RE()->charges_of_types);
        }},
       {"exclusion_range",
        [this](Variant const &v) {
          m_exclusion->do_set_parameter("exclusion_range", v);
        },
        [this]() { return m_exclusion->get_parameter("exclusion_range"); }},
       {"exclusion_radius_per_type",
        [this](Variant const &v) {
          m_exclusion->do_set_parameter("exclusion_radius_per_type", v);
        },
        [this]() {
          return m_exclusion->get_parameter("exclusion_radius_per_type");
        }}});
}

Variant ReactionAlgorithm::do_call_method(std::string const &name,
                                          VariantMap const &params) {
  if (name == "calculate_factorial_expression") {
    if (context()->is_head_node()) {
      auto &bookkeeping = RE()->get_old_system_state();
      auto &old_particle_numbers = bookkeeping.old_particle_numbers;
      auto &reaction = *m_reactions[bookkeeping.reaction_id]->get_reaction();
      return calculate_factorial_expression(reaction, old_particle_numbers);
    }
    return {};
  }
  if (name == "get_random_reaction_index") {
    return RE()->i_random(static_cast<int>(RE()->reactions.size()));
  }
  if (name == "potential_energy") {
    return RE()->calculate_potential_energy();
  }
  if (name == "create_new_trial_state") {
    auto const reaction_id = get_value<int>(params, "reaction_id");
    Variant result{};
    context()->parallel_try_catch([&]() {
      auto const optional = RE()->create_new_trial_state(reaction_id);
      if (optional) {
        result = *optional;
      }
    });
    return result;
  }
  if (name == "make_reaction_mc_move_attempt") {
    auto const bf = get_value<double>(params, "bf");
    auto const E_pot_old = get_value<double>(params, "E_pot_old");
    auto const E_pot_new = get_value<double>(params, "E_pot_new");
    auto const reaction_id = get_value<int>(params, "reaction_id");
    Variant result;
    context()->parallel_try_catch([&]() {
      result = RE()->make_reaction_mc_move_attempt(reaction_id, bf, E_pot_old,
                                                   E_pot_new);
    });
    return result;
  }
  if (name == "setup_bookkeeping_of_empty_pids") {
    RE()->setup_bookkeeping_of_empty_pids();
  } else if (name == "remove_constraint") {
    RE()->remove_constraint();
  } else if (name == "set_cylindrical_constraint_in_z_direction") {
    context()->parallel_try_catch([&]() {
      RE()->set_cyl_constraint(get_value<double>(params, "center_x"),
                               get_value<double>(params, "center_y"),
                               get_value<double>(params, "radius"));
    });
  } else if (name == "set_wall_constraints_in_z_direction") {
    context()->parallel_try_catch([&]() {
      RE()->set_slab_constraint(get_value<double>(params, "slab_start_z"),
                                get_value<double>(params, "slab_end_z"));
    });
  } else if (name == "get_wall_constraints_in_z_direction") {
    if (context()->is_head_node()) {
      return RE()->get_slab_constraint_parameters();
    }
    return {};
  } else if (name == "set_volume") {
    context()->parallel_try_catch(
        [&]() { RE()->set_volume(get_value<double>(params, "volume")); });
  } else if (name == "get_volume") {
    return RE()->get_volume();
  } else if (name == "get_acceptance_rate_reaction") {
    auto const index = get_value<int>(params, "reaction_id");
    context()->parallel_try_catch([&]() {
      if (index < 0 or index >= static_cast<int>(m_reactions.size())) {
        throw std::out_of_range("No reaction with id " + std::to_string(index));
      }
    });
    return m_reactions[index]->get_reaction()->get_acceptance_rate();
  } else if (name == "set_non_interacting_type") {
    auto const type = get_value<int>(params, "type");
    context()->parallel_try_catch([&]() {
      if (type < 0) {
        throw std::domain_error("Invalid type: " + std::to_string(type));
      }
    });
    RE()->non_interacting_type = type;
    ::init_type_map(type);
  } else if (name == "get_non_interacting_type") {
    return RE()->non_interacting_type;
  } else if (name == "displacement_mc_move_for_particles_of_type") {
    auto const type = get_value<int>(params, "type_mc");
    auto const n_particles =
        get_value_or<int>(params, "particle_number_to_be_changed", 1);
    auto result = false;
    context()->parallel_try_catch([&]() {
      result = RE()->make_displacement_mc_move_attempt(type, n_particles);
    });
    return result;
  } else if (name == "delete_particle") {
    context()->parallel_try_catch(
        [&]() { RE()->delete_particle(get_value<int>(params, "p_id")); });
  } else if (name == "delete_reaction") {
    auto const reaction_id = get_value<int>(params, "reaction_id");
    auto const index = get_reaction_index(reaction_id);
    // delete forward and backward reactions
    delete_reaction(index + 1);
    delete_reaction(index + 0);
  } else if (name == "add_reaction") {
    auto const reaction =
        get_value<std::shared_ptr<SingleReaction>>(params, "reaction");
    m_reactions.push_back(reaction);
    RE()->add_reaction(reaction->get_reaction());
  } else if (name == "change_reaction_constant") {
    auto const gamma = get_value<double>(params, "gamma");
    auto const reaction_id = get_value<int>(params, "reaction_id");
    context()->parallel_try_catch([&]() {
      if (reaction_id % 2 == 1) {
        throw std::invalid_argument("Only forward reactions can be selected");
      }
      if (gamma <= 0.) {
        throw std::domain_error("gamma needs to be a strictly positive value");
      }
    });
    auto const index = get_reaction_index(reaction_id);
    m_reactions[index + 0]->get_reaction()->gamma = gamma;
    m_reactions[index + 1]->get_reaction()->gamma = 1. / gamma;
  } else if (name == "set_charge_of_type") {
    auto const type = get_value<int>(params, "type");
    auto const charge = get_value<double>(params, "charge");
    RE()->charges_of_types[type] = charge;
  } else if (context()->is_head_node()) {
    throw std::runtime_error("unknown method '" + name + "()'");
  }
  return {};
}

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */
