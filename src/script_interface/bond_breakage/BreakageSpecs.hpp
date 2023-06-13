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

#pragma once

#include "BreakageSpec.hpp"

#include "core/bond_breakage/bond_breakage.hpp"

#include "script_interface/ObjectMap.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace ScriptInterface {
namespace BondBreakage {
class BreakageSpecs : public ObjectMap<BreakageSpec> {
  using container_type = std::unordered_map<int, std::shared_ptr<BreakageSpec>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

private:
  key_type insert_in_core(mapped_type const &) override {
    if (context()->is_head_node()) {
      throw std::runtime_error(
          "Inserting breakage spec without a bond type is not permitted.");
    }
    return {};
  }
  void insert_in_core(key_type const &key,
                      mapped_type const &obj_ptr) override {
    auto core_spec = obj_ptr->breakage_spec();
    ::BondBreakage::insert_spec(key, core_spec);
  }
  void erase_in_core(key_type const &key) override {
    ::BondBreakage::erase_spec(key);
  }
};
} // namespace BondBreakage
} // namespace ScriptInterface
