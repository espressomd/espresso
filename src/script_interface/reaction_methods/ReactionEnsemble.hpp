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

#ifndef SCRIPT_INTERFACE_REACTION_METHODS_REACTION_ENSEMBLE_HPP
#define SCRIPT_INTERFACE_REACTION_METHODS_REACTION_ENSEMBLE_HPP

#include "ReactionAlgorithm.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/reaction_methods/ReactionAlgorithm.hpp"
#include "core/reaction_methods/ReactionEnsemble.hpp"

#include <memory>

namespace ScriptInterface {
namespace ReactionMethods {

class ReactionEnsemble : public ReactionAlgorithm {
public:
  std::shared_ptr<::ReactionMethods::ReactionAlgorithm> RE() override {
    return m_re;
  }

  void do_construct(VariantMap const &params) override {
    m_re = std::make_shared<::ReactionMethods::ReactionEnsemble>(
        get_value<int>(params, "seed"), get_value<double>(params, "kT"),
        get_value<double>(params, "exclusion_radius"));
  }

private:
  std::shared_ptr<::ReactionMethods::ReactionEnsemble> m_re;
};
} /* namespace ReactionMethods */
} /* namespace ScriptInterface */

#endif