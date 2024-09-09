/*
 * Copyright (C) 2021-2024 The ESPResSo project
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

#include <config/config.hpp>

#ifdef COLLISION_DETECTION

#include "script_interface/interactions/BondedInteractions.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/collision_detection/ActiveProtocol.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>

namespace ScriptInterface::CollisionDetection {

class Protocol : public AutoParameters<Protocol> {
  std::weak_ptr<BondedInteractionsMap> m_bonded_ias;
  std::weak_ptr<Interactions::BondedInteractions> m_so_bonded_ias;
  std::unique_ptr<VariantMap> m_params;

public:
  virtual std::shared_ptr<::CollisionDetection::ActiveProtocol> protocol() = 0;
  void bind_bonded_ias(
      std::weak_ptr<BondedInteractionsMap> bonded_ias,
      std::weak_ptr<Interactions::BondedInteractions> so_bonded_ias) {
    m_bonded_ias = std::move(bonded_ias);
    m_so_bonded_ias = std::move(so_bonded_ias);
    if (m_params) {
      do_initialize(*m_params);
      m_params.reset();
    }
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

protected:
  auto find_bond_id(Variant const &v) const {
    auto bonded_ias = m_bonded_ias.lock();
    assert(bonded_ias && "This protocol is not bound to a system");
    std::optional<int> retval = std::nullopt;
    if (is_type<int>(v)) {
      auto const bond_id = get_value<int>(v);
      if (bonded_ias->contains(bond_id)) {
        retval = bond_id;
      }
    } else {
      auto obj = get_value<std::shared_ptr<Interactions::BondedInteraction>>(v);
      retval = bonded_ias->find_bond_id(obj->bonded_ia());
    }
    if (not retval) {
      throw std::invalid_argument("Bond in parameter list was "
                                  "not added to the system");
    }
    return *retval;
  }
  auto get_bond_variant_by_id(int bond_id) {
    auto so_bonded_ias = m_so_bonded_ias.lock();
    assert(so_bonded_ias != nullptr);
    return so_bonded_ias->do_call_method("get", {{"key", bond_id}});
  }
  virtual void do_initialize(VariantMap const &params) = 0;
};

} // namespace ScriptInterface::CollisionDetection

#endif // COLLISION_DETECTION
