/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

/**
 * @file
 * Reusable particle packing/unpacking for ghost communications.
 */

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "ParticleList.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace GhostComm {

namespace detail {
/**
 * @brief Allocator wrapper that skips value-initialization on grow.
 *
 * `std::vector<T, DefaultInitAllocator<T>>::resize(n)` extends storage
 * without zero-filling when n > size(), leaving new bytes in an
 * indeterminate state.  Use only when all grown bytes are overwritten
 * before being read (e.g. by pack_cells / MPI irecv).
 *
 * All other allocator operations (allocate, deallocate, copy/move) are
 * forwarded to `std::allocator<T>`.
 */
template <class T> struct DefaultInitAllocator : std::allocator<T> {
  using value_type = T;
  // Inherit all constructors from std::allocator<T>.
  using std::allocator<T>::allocator;

  template <class U> struct rebind {
    using other = DefaultInitAllocator<U>;
  };

  /**
   * @brief Called by vector::resize for each newly appended element.
   * Plain `new (p) T` uses default-initialization: trivial types like
   * `char` are left uninitialized rather than zero-filled.
   */
  void construct(T *p) noexcept { ::new (static_cast<void *>(p)) T; }

  /** Forward non-default construction unchanged. */
  template <class... Args> void construct(T *p, Args &&...args) {
    ::new (static_cast<void *>(p)) T(std::forward<Args>(args)...);
  }
};
} // namespace detail

/**
 * @brief Class that stores marshalled data for ghost communications.
 * To store and retrieve data, use the adapter classes below.
 *
 * Invariant: every byte in the non-bond buffer is written by the caller
 * (pack_cells / MPI irecv) before it is read (unpack_cells / add_forces).
 * The buffer therefore uses a default-init allocator so that resize()-on-grow
 * does not zero-fill bytes that will be immediately overwritten, eliminating
 * the `__memset_avx2` overhead visible in perf profiles at large particle
 * counts.  The bond buffer is kept as a plain std::vector<char> because it is
 * serialized via boost.mpi (which requires exact byte semantics) and is on the
 * cold resort-only path.
 */
class CommBuf {
  using ByteVec = std::vector<char, detail::DefaultInitAllocator<char>>;

public:
  /** Returns a pointer to the non-bond storage. */
  char *data() { return buf.data(); }
  const char *data() const { return buf.data(); }

  /** Returns the number of elements in the non-bond storage. */
  std::size_t size() const { return buf.size(); }

  /**
   * @brief Resizes the underlying storage s.t. the object is capable of holding
   * @p new_size chars.  Bytes added by growth are left uninitialized;
   * the caller must write them (pack_cells / MPI irecv) before reading.
   */
  void resize(std::size_t new_size) { buf.resize(new_size); }

  /** Returns a reference to the bond storage. */
  auto &bonds() { return bondbuf; }
  const auto &bonds() const { return bondbuf; }

  auto make_span() { return std::span(buf.data(), buf.size()); }

private:
  ByteVec buf; ///< Buffer for everything but bonds (default-init grow)
  std::vector<char>
      bondbuf; ///< Buffer for bond lists (boost.mpi-serialized, cold path)
};

/**
 * @brief Calculate the per-particle transmit size for the given data parts.
 */
std::size_t calc_transmit_size(BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Calculate the total transmit size for a list of cells and data parts.
 *
 * When GHOSTTRANS_PARTNUM is set, returns sizeof(unsigned int) per cell.
 * Otherwise returns the total number of particles times the per-particle size.
 */
std::size_t calc_transmit_size(std::span<ParticleList *const> cells,
                               BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Pack particle data from cells into a communication buffer.
 *
 * Handles GHOSTTRANS_PARTNUM (writes cell sizes) and GHOSTTRANS_BONDS
 * (separate bond buffer) special cases.
 *
 * @param buf        Buffer to pack into (resized as needed).
 * @param cells      Source particle lists.
 * @param shift      Ghost shift applied to particle positions on save.
 * @param box_geo    Box geometry for fold_position.
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void pack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                unsigned data_parts);

/**
 * @brief Unpack particle data from a communication buffer into cells.
 *
 * Handles GHOSTTRANS_PARTNUM (resizes ghost cells) and GHOSTTRANS_BONDS
 * (separate bond buffer) special cases.
 *
 * @param buf        Buffer to unpack from.
 * @param cells      Destination particle lists.
 * @param box_geo    Box geometry (unused here, kept for API symmetry).
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void unpack_cells(CommBuf &buf, std::span<ParticleList *const> cells,
                  BoxGeometry const &box_geo, unsigned data_parts);

/**
 * @brief Add forces (and optionally torques) from a communication buffer
 *        to particles in cells.
 *
 * @param buf        Buffer produced by a force-reduce exchange.
 * @param cells      Destination (owned) particle lists.
 * @param data_parts Bitmask of GHOSTTRANS_* flags; must include
 *                   GHOSTTRANS_FORCE, may include GHOSTTRANS_TORQUE.
 */
void add_forces(CommBuf &buf, std::span<ParticleList *const> cells,
                unsigned data_parts);

#ifdef ESPRESSO_BOND_CONSTRAINT
/**
 * @brief Add rattle corrections from a communication buffer to particles.
 */
void add_rattle(CommBuf &buf, std::span<ParticleList *const> cells);
#endif

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
/**
 * @brief Add dipole fields from a communication buffer to particles.
 */
void add_dip_fld(CommBuf &buf, std::span<ParticleList *const> cells);
#endif

/**
 * @brief Copy particle data from src to dst applying a ghost shift.
 *
 * For GHOSTTRANS_PARTNUM, resizes dst to match src.
 * Otherwise serializes each particle from src and deserializes into dst,
 * applying the shift and folding the position. Bond data is copied directly.
 *
 * @param src        Source particle list.
 * @param dst        Destination particle list.
 * @param shift      Ghost shift to apply.
 * @param box_geo    Box geometry for fold_position.
 * @param data_parts Bitmask of GHOSTTRANS_* flags.
 */
void local_cell_copy(ParticleList &src, ParticleList &dst,
                     Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                     unsigned data_parts);

} // namespace GhostComm
