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

#include "initialize.hpp"
#include "BondBreakage.hpp"
#include "BreakageSpec.hpp"
#include "BreakageSpecs.hpp"

namespace ScriptInterface {
namespace BondBreakage {
void initialize(Utils::Factory<ObjectHandle> *f) {
  f->register_new<BreakageSpec>("BondBreakage::BreakageSpec");
  f->register_new<BreakageSpecs>("BondBreakage::BreakageSpecs");
  f->register_new<BondBreakage>("BondBreakage::BondBreakage");
}
} // namespace BondBreakage
} // namespace ScriptInterface
