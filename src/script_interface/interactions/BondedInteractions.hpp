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

#pragma once

#include "BondedInteraction.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"

#include "script_interface/ObjectMap.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/system/Leaf.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Interactions {

using BondedInteractionsBase_t = ObjectMap<
    BondedInteraction,
    AutoParameters<ObjectMap<BondedInteraction, System::Leaf>, System::Leaf>>;

class BondedInteractions : public BondedInteractionsBase_t {
  using Base = BondedInteractionsBase_t;
  using container_type =
      std::unordered_map<int, std::shared_ptr<BondedInteraction>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

private:
  container_type m_bonds;
  std::shared_ptr<::BondedInteractionsMap> m_handle;
  std::unique_ptr<VariantMap> m_params;

  void do_construct(VariantMap const &params) override {
    m_handle = std::make_shared<::BondedInteractionsMap>();
    m_params = std::make_unique<VariantMap>(params);
  }

  void on_bind_system(::System::System &system) override {
    system.bonded_ias = m_handle;
    m_handle->bind_system(m_system.lock());
    m_handle->on_ia_change();
    if (m_params and not m_params->empty()) {
      restore_from_checkpoint(*m_params);
    }
    m_params.reset();
  }

  key_type insert_in_core(mapped_type const &obj_ptr) override {
    key_type key{};
    context()->parallel_try_catch(
        [&]() { key = m_handle->insert(obj_ptr->bonded_ia()); });
    m_bonds[key] = std::move(obj_ptr);
    return key;
  }

  void insert_in_core(key_type const &key,
                      mapped_type const &obj_ptr) override {
    context()->parallel_try_catch(
        [&]() { m_handle->insert(key, obj_ptr->bonded_ia()); });
    m_bonds[key] = std::move(obj_ptr);
  }

  void erase_in_core(key_type const &key) override {
    m_handle->erase(key);
    m_bonds.erase(key);
  }

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_size") {
      return {static_cast<int>(m_handle->size())};
    }

    if (name == "get_bond_ids") {
      std::vector<int> bond_ids;
      for (auto const &kv : *m_handle)
        bond_ids.push_back(kv.first);
      return bond_ids;
    }

    if (name == "has_bond") {
      auto const bond_id = get_key(params.at("bond_id"));
      return {m_bonds.contains(bond_id)};
    }

    if (name == "get_bond") {
      auto const bond_id = get_key(params.at("bond_id"));
      // core and script interface must agree
      assert(m_bonds.count(bond_id) == m_handle->count(bond_id));
      if (not context()->is_head_node())
        return {};
      // bond must exist
      if (not m_bonds.contains(bond_id)) {
        throw std::out_of_range("The bond with id " + std::to_string(bond_id) +
                                " is not yet defined.");
      }
      return {m_bonds.at(bond_id)};
    }

    if (name == "get_zero_based_type") {
      auto const bond_id = get_key(params.at("bond_id"));
      return m_handle->get_zero_based_type(bond_id);
    }

    return Base::do_call_method(name, params);
  }
};

} // namespace Interactions
} // namespace ScriptInterface
