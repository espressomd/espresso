/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "LatticeIndices.hpp"

#include <script_interface/ScriptInterface.hpp>

#include <walberla_bridge/LatticeWalberla.hpp>

#include <utils/Vector.hpp>

#include <cassert>
#include <tuple>
#include <vector>

namespace ScriptInterface::walberla {

template <class FieldSerializer> class LatticeSlice : public LatticeIndices {
protected:
  Utils::Vector3i m_slice_lower_corner;
  Utils::Vector3i m_slice_upper_corner;
  std::vector<int> m_shape;

public:
  virtual ::LatticeWalberla const &get_lattice() const = 0;

private:
  auto get_sentinel_index(::LatticeWalberla const &lattice) const {
    return -(static_cast<int>(lattice.get_ghost_layers()) + 1);
  }

  auto get_slices_bounding_boxes() const {
    auto const &lattice = get_lattice();
    auto const &slice_lower_corner = m_slice_lower_corner;
    auto const &slice_upper_corner = m_slice_upper_corner;
    assert(slice_upper_corner <= lattice.get_grid_dimensions());
    assert(slice_lower_corner >= Utils::Vector3i::broadcast(0));
    auto const sentinel = get_sentinel_index(lattice);
    auto [local_lower_corner, local_upper_corner] =
        lattice.get_local_grid_range();
    for (auto const i : {0, 1, 2}) {
      if (local_lower_corner[i] >= slice_upper_corner[i] or
          slice_lower_corner[i] >= local_upper_corner[i]) {
        local_lower_corner[i] = sentinel;
        local_upper_corner[i] = sentinel;
      } else {
        if (slice_lower_corner[i] > local_lower_corner[i]) {
          local_lower_corner[i] = slice_lower_corner[i];
        }
        if (slice_upper_corner[i] < local_upper_corner[i]) {
          local_upper_corner[i] = slice_upper_corner[i];
        }
      }
    }
    return std::make_tuple(slice_lower_corner, slice_upper_corner,
                           local_lower_corner, local_upper_corner);
  }

protected:
  template <class LatticeModel, typename T>
  Variant gather_3d(VariantMap const &params, std::vector<int> const &data_dims,
                    LatticeModel const &lattice_model,
                    std::vector<T> (LatticeModel::*getter)(
                        Utils::Vector3i const &, Utils::Vector3i const &) const,
                    double units_conversion = 1.) const;

  template <class LatticeModel, typename T>
  void scatter_3d(VariantMap const &params, std::vector<int> const &data_dims,
                  LatticeModel &lattice_model,
                  void (LatticeModel::*setter)(Utils::Vector3i const &,
                                               Utils::Vector3i const &,
                                               std::vector<T> const &),
                  double units_conversion = 1.);
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
