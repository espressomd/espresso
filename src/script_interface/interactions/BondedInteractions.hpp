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

#ifndef SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTIONS_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTIONS_HPP

#include "BondedInteraction.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"

#include "script_interface/ObjectMap.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Interactions {
class BondedInteractions : public ObjectMap<BondedInteraction> {
  using container_type =
      std::unordered_map<int, std::shared_ptr<BondedInteraction>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;

  BondedInteractions() : ObjectMap<BondedInteraction>::ObjectMap() {
    add_parameters({
        {"_objects", AutoParameter::read_only,
         []() {
           // deactivate serialization (done at the Python level)
           return make_unordered_map_of_variants(container_type{});
         }},
    });
  }

private:
  container_type m_bonds;

  key_type insert_in_core(mapped_type const &obj_ptr) override {
    auto const key = ::bonded_ia_params.insert(obj_ptr->bonded_ia());
    m_bonds[key] = std::move(obj_ptr);
    mpi_update_cell_system_ia_range_local();
    return key;
  }

  void insert_in_core(key_type const &key,
                      mapped_type const &obj_ptr) override {
    ::bonded_ia_params.insert(key, obj_ptr->bonded_ia());
    m_bonds[key] = std::move(obj_ptr);
    mpi_update_cell_system_ia_range_local();
  }

  void erase_in_core(key_type const &key) override {
    ::bonded_ia_params.erase(key);
    m_bonds.erase(key);
    mpi_update_cell_system_ia_range_local();
  }

protected:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "get_size") {
      return {static_cast<int>(::bonded_ia_params.size())};
    }

    if (name == "get_bond_ids") {
      std::vector<int> bond_ids;
      for (auto const &kv : ::bonded_ia_params)
        bond_ids.push_back(kv.first);
      return bond_ids;
    }

    if (name == "has_bond") {
      auto const bond_id = get_key(params.at("bond_id"));
      return {m_bonds.count(bond_id) != 0};
    }

    if (name == "get_bond") {
      auto const bond_id = get_key(params.at("bond_id"));
      // core and script interface must agree
      assert(m_bonds.count(bond_id) == ::bonded_ia_params.count(bond_id));
      if (not context()->is_head_node())
        return {};
      // bond must exist
      if (m_bonds.count(bond_id) == 0) {
        throw std::out_of_range("The bond with id " + std::to_string(bond_id) +
                                " is not yet defined.");
      }
      return {m_bonds.at(bond_id)};
    }

    if (name == "get_zero_based_type") {
      auto const bond_id = get_key(params.at("bond_id"));
      return ::bonded_ia_params.get_zero_based_type(bond_id);
    }

    return ObjectMap<BondedInteraction>::do_call_method(name, params);
  }
};
} // namespace Interactions
} // namespace ScriptInterface

#endif
