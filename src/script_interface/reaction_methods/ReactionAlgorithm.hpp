/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include "ExclusionRadius.hpp"
#include "SingleReaction.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/particle_data/ParticleHandle.hpp"

#include "core/reaction_methods/ReactionAlgorithm.hpp"
#include "core/reaction_methods/utils.hpp"

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
  std::shared_ptr<ExclusionRadius> m_exclusion;
  /**
   * Check reaction id is within the reaction container bounds.
   * Since each reaction has a corresponding backward reaction,
   * the total number of reactions is doubled. Return the
   * corresponding index for @ref ReactionAlgorithm::m_reactions.
   */
  int get_reaction_index(int reaction_id) const {
    auto const index = 2 * reaction_id;
    context()->parallel_try_catch([&]() {
      if (index < 0 or index >= static_cast<int>(m_reactions.size())) {
        throw std::out_of_range("No reaction with id " +
                                std::to_string(reaction_id));
      }
    });
    return index;
  }

  void setup_neighbor_search(VariantMap const &params) {
    auto so_ptr = get_value<ObjectRef>(params, "exclusion");
    m_exclusion = std::dynamic_pointer_cast<ExclusionRadius>(so_ptr);
    m_exclusion->do_set_parameter("search_algorithm",
                                  Variant{get_value_or<std::string>(
                                      params, "search_algorithm", "order_n")});
  }

public:
  virtual std::shared_ptr<::ReactionMethods::ReactionAlgorithm> RE() = 0;
  virtual std::shared_ptr<::ReactionMethods::ReactionAlgorithm> const
  RE() const = 0;

  ReactionAlgorithm();

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

protected:
  virtual double calculate_factorial_expression(
      ::ReactionMethods::SingleReaction const &reaction,
      std::unordered_map<int, int> const &particle_numbers) const {
    return ::ReactionMethods::calculate_factorial_expression(reaction,
                                                             particle_numbers);
  }

private:
  void delete_reaction(int reaction_id) {
    m_reactions.erase(m_reactions.begin() + reaction_id);
    RE()->delete_reaction(reaction_id);
  }

  std::string get_internal_state() const override {
    assert(not RE()->is_reaction_under_way());
    throw std::runtime_error("Reaction methods do not support checkpointing");
  }
};

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif
