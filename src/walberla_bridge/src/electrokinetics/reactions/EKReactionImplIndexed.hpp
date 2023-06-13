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

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactant.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBase.hpp>

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace walberla {

class EKReactionImplIndexed : public EKReactionBase {
private:
  std::size_t m_flagfield_id;
  std::size_t m_indexvector_id;

  bool m_pending_changes;

public:
  EKReactionImplIndexed(std::shared_ptr<LatticeWalberla> lattice,
                        std::vector<std::shared_ptr<EKReactant>> reactants,
                        double coefficient);
  ~EKReactionImplIndexed() override = default;

  using EKReactionBase::get_coefficient;
  using EKReactionBase::get_lattice;
  using EKReactionBase::get_reactants;

  void perform_reaction() override;

  void set_node_is_boundary(Utils::Vector3i const &node, bool is_boundary);
  [[nodiscard]] boost::optional<bool>
  get_node_is_boundary(Utils::Vector3i const &node);

  [[nodiscard]] auto get_indexvector_id() const noexcept {
    return m_indexvector_id;
  }
  [[nodiscard]] auto get_flagfield_id() const noexcept {
    return m_flagfield_id;
  }

  void boundary_update();
};

} // namespace walberla
