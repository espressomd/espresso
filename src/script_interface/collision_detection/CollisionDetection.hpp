/*
 * Copyright (C) 2010-2024 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include "Protocol.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/interactions/BondedInteraction.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"
#include "script_interface/system/Leaf.hpp"
#include "script_interface/system/System.hpp"

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "core/collision_detection/CollisionDetection.hpp"

#include <memory>

namespace ScriptInterface::CollisionDetection {

class CollisionDetection
    : public AutoParameters<CollisionDetection, System::Leaf> {
  std::shared_ptr<::CollisionDetection::CollisionDetection> m_handle;
  std::unique_ptr<VariantMap> m_params;
  std::shared_ptr<Protocol> m_protocol;
  std::weak_ptr<::BondedInteractionsMap> m_bonded_ias;
  std::weak_ptr<Interactions::BondedInteractions> m_so_bonded_ias;

public:
  CollisionDetection() {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            if (is_none(value)) {
              m_protocol = nullptr;
              m_handle->unset_protocol();
              return;
            }
            auto const m_protocol_backup = m_protocol;
            try {
              context()->parallel_try_catch([&]() {
                m_protocol = get_value<std::shared_ptr<Protocol>>(value);
                m_protocol->bind_bonded_ias(m_bonded_ias, m_so_bonded_ias);
                m_handle->set_protocol(m_protocol->protocol());
              });
            } catch (...) {
              m_protocol = m_protocol_backup;
              if (m_protocol) {
                m_handle->set_protocol(m_protocol->protocol());
              } else {
                m_protocol = nullptr;
                m_handle->unset_protocol();
              }
              throw;
            }
          },
          [this]() {
            if (m_protocol)
              return make_variant(m_protocol);
            return make_variant(none);
          }}});
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

  void on_bind_system(::System::System &system) override {
    m_handle = system.collision_detection;
    m_handle->bind_system(m_system.lock());
    m_bonded_ias = system.bonded_ias;
  }

  void attach(std::weak_ptr<Interactions::BondedInteractions> bonded_ias) {
    m_so_bonded_ias = bonded_ias;
    if (m_params) {
      AutoParameters<CollisionDetection, System::Leaf>::do_construct(*m_params);
      m_params.reset();
    }
  }
};

} // namespace ScriptInterface::CollisionDetection

#endif // COLLISION_DETECTION
