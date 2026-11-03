/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#ifdef ESPRESSO_WALBERLA

#include "EKReaction.hpp"
#include "EKSpeciesSlice.hpp"
#include "LatticeSlice.hpp"

#include <script_interface/ScriptInterface.hpp>

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/electrokinetics/reactions/EKReactionBaseIndexed.hpp>

#include <utils/Vector.hpp>

#include <cassert>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {

/**
 * @brief Adapter that exposes @c EKReactionBaseIndexed with the
 *        interface expected by @c LatticeSlice::gather_3d / scatter_3d.
 *
 * @c LatticeSlice requires the lattice model to provide
 * @c get_lattice() returning a @c LatticeWalberla const &, but
 * @c EKReactionBase::get_lattice() returns a @c shared_ptr.
 * This adapter bridges that mismatch.
 */
struct EKReactionSliceAdapter {
  ::walberla::EKReactionBaseIndexed *m_impl;

  [[nodiscard]] ::LatticeWalberla const &get_lattice() const noexcept {
    return *m_impl->get_lattice();
  }

  [[nodiscard]] std::vector<int>
  get_slice_is_boundary(Utils::Vector3i const &lower_corner,
                        Utils::Vector3i const &upper_corner) const {
    return m_impl->get_slice_is_boundary(lower_corner, upper_corner);
  }

  void set_slice_is_boundary(Utils::Vector3i const &lower_corner,
                             Utils::Vector3i const &upper_corner,
                             std::vector<int> const &values) {
    m_impl->set_slice_is_boundary(lower_corner, upper_corner, values);
  }

  void ghost_communication() { m_impl->ghost_communication(); }
};

class EKIndexedReactionSlice : public LatticeSlice<EKFieldSerializer> {
  std::shared_ptr<::walberla::EKReactionBaseIndexed> m_reaction_impl;
  std::shared_ptr<EKIndexedReaction> m_reaction_sip;
  std::unordered_map<std::string, std::vector<int>> m_shape_val;

public:
  void do_construct(VariantMap const &params) override {
    m_reaction_sip =
        get_value<std::shared_ptr<EKIndexedReaction>>(params, "parent_sip");
    m_reaction_impl = m_reaction_sip->get_impl();
    assert(m_reaction_impl);
    m_shape = get_value<std::vector<int>>(params, "shape");
    m_slice_lower_corner =
        get_value<Utils::Vector3i>(params, "slice_lower_corner");
    m_slice_upper_corner =
        get_value<Utils::Vector3i>(params, "slice_upper_corner");
    m_shape_val["is_boundary"] = std::vector<int>(1, 1);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  ::LatticeWalberla const &get_lattice() const override {
    return *m_reaction_impl->get_lattice();
  }
};

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
