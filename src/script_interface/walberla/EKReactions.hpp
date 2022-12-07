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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTIONS_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_EK_REACTIONS_HPP

#include "config/config.hpp"

#ifdef WALBERLA

#include "EKReaction.hpp"

#include "core/grid_based_algorithms/ek_reactions.hpp"

#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface::walberla {

class EKReactions : public ObjectList<EKReaction> {
  void add_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {
    EK::ek_reactions.add(obj_ptr->get_instance());
  }
  void remove_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {
    EK::ek_reactions.remove(obj_ptr->get_instance());
  }
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
