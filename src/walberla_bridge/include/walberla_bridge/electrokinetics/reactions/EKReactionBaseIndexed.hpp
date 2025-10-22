/*
 * Copyright (C) 2024 The ESPResSo project
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

#include "EKReactionBase.hpp"

#include <utils/Vector.hpp>

#include <optional>

namespace walberla {

class EKReactionBaseIndexed : public EKReactionBase {
public:
  using EKReactionBase::EKReactionBase;
  ~EKReactionBaseIndexed() override = default;
  virtual void set_node_is_boundary(Utils::Vector3i const &node,
                                    bool is_boundary) = 0;
  virtual std::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node) = 0;
};

} // namespace walberla
