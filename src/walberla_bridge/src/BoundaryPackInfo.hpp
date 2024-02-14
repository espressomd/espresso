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

#include <core/debug/Debug.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>
#include <domain_decomposition/IBlock.h>
#include <field/communication/PackInfo.h>
#include <stencil/Directions.h>

#include <memory>
#include <tuple>
#include <utility>

namespace walberla {
namespace field {
namespace communication {

template <typename GhostLayerField_T, typename Boundary_T>
class BoundaryPackInfo : public PackInfo<GhostLayerField_T> {
protected:
  using PackInfo<GhostLayerField_T>::bdId_;

public:
  using PackInfo<GhostLayerField_T>::PackInfo;
  using PackInfo<GhostLayerField_T>::numberOfGhostLayersToCommunicate;

  ~BoundaryPackInfo() override = default;

  void setup_boundary_handle(std::shared_ptr<LatticeWalberla> lattice,
                             std::shared_ptr<Boundary_T> boundary) {
    m_lattice = std::move(lattice);
    m_boundary = std::move(boundary);
  }

  bool constantDataExchange() const override { return false; }
  bool threadsafeReceiving() const override { return true; }

  void communicateLocal(IBlock const *sender, IBlock *receiver,
                        stencil::Direction dir) override {
    mpi::SendBuffer sBuffer;
    packDataImpl(sender, dir, sBuffer);
    mpi::RecvBuffer rBuffer(sBuffer);
    unpackData(receiver, stencil::inverseDir[dir], rBuffer);
  }

  void unpackData(IBlock *receiver, stencil::Direction dir,
                  mpi::RecvBuffer &buffer) override {

    auto *flag_field = receiver->getData<GhostLayerField_T>(bdId_);
    WALBERLA_ASSERT_NOT_NULLPTR(flag_field);
    WALBERLA_ASSERT_NOT_NULLPTR(m_boundary);
    WALBERLA_ASSERT_NOT_NULLPTR(m_lattice);

    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    auto const gl = numberOfGhostLayersToCommunicate(flag_field);
    auto const begin = [gl, dir](auto const *flag_field) {
      return flag_field->beginGhostLayerOnly(gl, dir);
    };

#ifndef NDEBUG
    uint_t xSize, ySize, zSize, bSize;
    buffer >> xSize >> ySize >> zSize >> bSize;
    uint_t buf_size{0u};
    for (auto it = begin(flag_field); it != flag_field->end(); ++it) {
      if (isFlagSet(it, boundary_flag)) {
        ++buf_size;
      }
    }
    WALBERLA_ASSERT_EQUAL(xSize, flag_field->xSize());
    WALBERLA_ASSERT_EQUAL(ySize, flag_field->ySize());
    WALBERLA_ASSERT_EQUAL(zSize, flag_field->zSize());
    WALBERLA_ASSERT_EQUAL(bSize, buf_size);
#endif

    auto const offset = std::get<0>(m_lattice->get_local_grid_range());
    typename Boundary_T::value_type value;
    for (auto it = begin(flag_field); it != flag_field->end(); ++it) {
      if (isFlagSet(it, boundary_flag)) {
        auto const node = offset + Utils::Vector3i{{it.x(), it.y(), it.z()}};
        buffer >> value;
        m_boundary->unpack_node(node, value);
      }
    }
  }

protected:
  void packDataImpl(IBlock const *sender, stencil::Direction dir,
                    mpi::SendBuffer &buffer) const override {

    auto const *flag_field = sender->getData<GhostLayerField_T>(bdId_);
    WALBERLA_ASSERT_NOT_NULLPTR(flag_field);
    WALBERLA_ASSERT_NOT_NULLPTR(m_boundary);
    WALBERLA_ASSERT_NOT_NULLPTR(m_lattice);

    auto const boundary_flag = flag_field->getFlag(Boundary_flag);
    auto const gl = numberOfGhostLayersToCommunicate(flag_field);
    auto const begin = [gl, dir](auto const *flag_field) {
      return flag_field->beginSliceBeforeGhostLayer(dir, gl);
    };

#ifndef NDEBUG
    uint_t buf_size{0u};
    for (auto it = begin(flag_field); it != flag_field->end(); ++it) {
      if (isFlagSet(it, boundary_flag)) {
        ++buf_size;
      }
    }
    buffer << flag_field->xSize() << flag_field->ySize() << flag_field->zSize()
           << buf_size;
#endif

    auto const offset = std::get<0>(m_lattice->get_local_grid_range());
    for (auto it = begin(flag_field); it != flag_field->end(); ++it) {
      if (isFlagSet(it, boundary_flag)) {
        auto const node = offset + Utils::Vector3i{{it.x(), it.y(), it.z()}};
        buffer << m_boundary->get_node_value_at_boundary(node);
      }
    }
  }

private:
  std::shared_ptr<LatticeWalberla> m_lattice;
  std::shared_ptr<Boundary_T> m_boundary;
};

} // namespace communication
} // namespace field
} // namespace walberla
