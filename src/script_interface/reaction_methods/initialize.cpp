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

#include "initialize.hpp"

#include "ExclusionRadius.hpp"
#include "SingleReaction.hpp"

#include "ConstantpHEnsemble.hpp"
#include "ReactionEnsemble.hpp"
#include "WidomInsertion.hpp"

#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface {
namespace ReactionMethods {
void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<SingleReaction>("ReactionMethods::SingleReaction");
  om->register_new<WidomInsertion>("ReactionMethods::WidomInsertion");
  om->register_new<ReactionEnsemble>("ReactionMethods::ReactionEnsemble");
  om->register_new<ConstantpHEnsemble>("ReactionMethods::ConstantpHEnsemble");
  om->register_new<ExclusionRadius>("ReactionMethods::ExclusionRadius");
}
} // namespace ReactionMethods
} // namespace ScriptInterface
