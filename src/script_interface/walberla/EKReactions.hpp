/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>

#include "core/ek/EKReactions.hpp"
#include "core/ek/EKWalberla.hpp"

#include "EKReaction.hpp"

#include <script_interface/ObjectList.hpp>
#include <script_interface/ScriptInterface.hpp>

#include <memory>

namespace ScriptInterface::walberla {

class EKReactions : public ObjectList<EKReaction> {
  std::shared_ptr<::EK::EKWalberla::ek_reactions_type> m_ek_reactions;

  void add_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {
    m_ek_reactions->add(obj_ptr->get_instance());
  }
  void remove_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {
    m_ek_reactions->remove(obj_ptr->get_instance());
  }

protected:
  void do_construct(VariantMap const &params) override {
    m_ek_reactions = std::make_shared<::EK::EKWalberla::ek_reactions_type>();
  }

public:
  auto &get_handle() { return m_ek_reactions; }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
