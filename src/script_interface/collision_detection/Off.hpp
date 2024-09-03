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

#include "Protocol.hpp"

#include "core/collision_detection/Off.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>

namespace ScriptInterface::CollisionDetection {

class Off : public Protocol {
public:
  Off()
      : m_protocol{std::make_shared<::CollisionDetection::ActiveProtocol>(
            ::CollisionDetection::Off())} {}
  std::shared_ptr<::CollisionDetection::ActiveProtocol> protocol() override {
    return m_protocol;
  }

protected:
  void do_initialize(VariantMap const &) override {}

private:
  std::shared_ptr<::CollisionDetection::ActiveProtocol> m_protocol;
};

} // namespace ScriptInterface::CollisionDetection
#endif // COLLISION_DETECTION
