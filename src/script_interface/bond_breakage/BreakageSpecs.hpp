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
#include "script_interface/system/Leaf.hpp"

#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace ScriptInterface {
namespace BondBreakage {

class BreakageSpecs
    : public ObjectMap<
          BreakageSpec,
          AutoParameters<ObjectMap<BreakageSpec, System::Leaf>, System::Leaf>> {
  using container_type = std::unordered_map<int, std::shared_ptr<BreakageSpec>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

private:
  std::shared_ptr<::BondBreakage::BondBreakage> m_bond_breakage;

  void do_construct(VariantMap const &params) override {
    m_bond_breakage = std::make_shared<::BondBreakage::BondBreakage>();
    restore_from_checkpoint(params);
  }

  void on_bind_system(::System::System &system) override {
    system.bond_breakage = m_bond_breakage;
  }

  key_type insert_in_core(mapped_type const &) override {
    if (context()->is_head_node()) {
      throw std::runtime_error(
          "Inserting breakage spec without a bond type is not permitted.");
    }
    return {};
  }
  void insert_in_core(key_type const &key,
                      mapped_type const &obj_ptr) override {
    m_bond_breakage->breakage_specs.emplace(key, obj_ptr->breakage_spec());
  }
  void erase_in_core(key_type const &key) override {
    m_bond_breakage->breakage_specs.erase(key);
  }
};
} // namespace BondBreakage
} // namespace ScriptInterface
