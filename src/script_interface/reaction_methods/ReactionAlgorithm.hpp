/*
 * Copyright (C) 2021 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_REACTION_METHODS_REACTION_ALGORITHM_HPP
#define SCRIPT_INTERFACE_REACTION_METHODS_REACTION_ALGORITHM_HPP

#include "SingleReaction.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/reaction_methods/ReactionAlgorithm.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface {
namespace ReactionMethods {

class ReactionAlgorithm : public AutoParameters<ReactionAlgorithm> {
protected:
  /** Keep track of the script interface pointer of each reaction. */
  std::vector<std::shared_ptr<SingleReaction>> m_reactions;
  /**
   * Check reaction id is within the reaction container bounds.
   * Since each reaction has a corresponding backward reaction,
   * the total number of reactions is doubled. Return the
   * corresponding index for @ref ReactionAlgorithm::m_reactions.
   */
  int get_reaction_index(int reaction_id) const {
    auto const index = 2 * reaction_id;
    if (index < 0 or index >= static_cast<int>(m_reactions.size())) {
      throw std::out_of_range("This reaction is not present");
    }
    return index;
  }

public:
  virtual std::shared_ptr<::ReactionMethods::ReactionAlgorithm> RE() = 0;

  ReactionAlgorithm() {
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
            auto const key = get_value<std::string>(v);
            if (key == "order_n") {
              RE()->neighbor_search_order_n = true;
            } else if (key == "parallel") {
              RE()->neighbor_search_order_n = false;
            } else {
              throw std::invalid_argument("Unknown search algorithm '" + key +
                                          "'");
            }
          },
          [this]() {
            if (RE()->neighbor_search_order_n) {
              return std::string("order_n");
            }
            return std::string("parallel");
          }},
         {"exclusion_range", AutoParameter::read_only,
          [this]() { return RE()->get_exclusion_range(); }},
         {"exclusion_radius_per_type",
          [this](Variant const &v) {
            RE()->set_exclusion_radius_per_type(
                get_value<std::unordered_map<int, double>>(v));
          },
          [this]() {
            return make_unordered_map_of_variants(
                RE()->exclusion_radius_per_type);
          }}});
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override {
    if (name == "remove_constraint") {
      RE()->remove_constraint();
    } else if (name == "set_cylindrical_constraint_in_z_direction") {
      RE()->set_cyl_constraint(get_value<double>(parameters, "center_x"),
                               get_value<double>(parameters, "center_y"),
                               get_value<double>(parameters, "radius"));
    } else if (name == "set_wall_constraints_in_z_direction") {
      RE()->set_slab_constraint(get_value<double>(parameters, "slab_start_z"),
                                get_value<double>(parameters, "slab_end_z"));
    } else if (name == "get_wall_constraints_in_z_direction") {
      return RE()->get_slab_constraint_parameters();
    } else if (name == "set_volume") {
      RE()->set_volume(get_value<double>(parameters, "volume"));
    } else if (name == "get_volume") {
      return RE()->get_volume();
    } else if (name == "get_acceptance_rate_reaction") {
      auto const index = get_value<int>(parameters, "reaction_id");
      if (index < 0 or index >= static_cast<int>(m_reactions.size())) {
        throw std::out_of_range("This reaction is not present");
      }
      return m_reactions[index]->get_reaction()->get_acceptance_rate();
    } else if (name == "set_non_interacting_type") {
      RE()->non_interacting_type = get_value<int>(parameters, "type");
    } else if (name == "get_non_interacting_type") {
      return RE()->non_interacting_type;
    } else if (name == "reaction") {
      RE()->do_reaction(get_value_or<int>(parameters, "reaction_steps", 1));
    } else if (name == "displacement_mc_move_for_particles_of_type") {
      return RE()->do_global_mc_move_for_particles_of_type(
          get_value<int>(parameters, "type_mc"),
          get_value_or<int>(parameters, "particle_number_to_be_changed", 1));
    } else if (name == "check_reaction_method") {
      RE()->check_reaction_method();
    } else if (name == "delete_particle") {
      RE()->delete_particle(get_value<int>(parameters, "p_id"));
    } else if (name == "delete_reaction") {
      auto const reaction_id = get_value<int>(parameters, "reaction_id");
      auto const index = get_reaction_index(reaction_id);
      // delete forward and backward reactions
      delete_reaction(index + 1);
      delete_reaction(index + 0);
    } else if (name == "add_reaction") {
      auto const reaction =
          get_value<std::shared_ptr<SingleReaction>>(parameters, "reaction");
      m_reactions.push_back(reaction);
      RE()->add_reaction(reaction->get_reaction());
    } else if (name == "change_reaction_constant") {
      auto const gamma = get_value<double>(parameters, "gamma");
      auto const reaction_id = get_value<int>(parameters, "reaction_id");
      auto const index = get_reaction_index(reaction_id);
      m_reactions[index + 0]->get_reaction()->gamma = gamma;
      m_reactions[index + 1]->get_reaction()->gamma = 1. / gamma;
    } else if (name == "set_charge_of_type") {
      auto const type = get_value<int>(parameters, "type");
      auto const charge = get_value<double>(parameters, "charge");
      RE()->charges_of_types[type] = charge;
    } else {
      throw std::runtime_error(("unknown method '" + name + "()'").c_str());
    }
    return {};
  };

private:
  void delete_reaction(int reaction_id) {
    m_reactions.erase(m_reactions.begin() + reaction_id);
    RE()->delete_reaction(reaction_id);
  }

  std::string get_internal_state() const override {
    throw std::runtime_error("Reaction methods do not support checkpointing");
  }
};
} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif