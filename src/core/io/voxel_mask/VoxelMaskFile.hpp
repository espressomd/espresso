/*
 * Copyright (C) 2026 The ESPResSo project
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

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/communicator.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace IO {

/**
 * @brief Parallel read access to a 3D voxel mask stored in an HDF5 dataset.
 *
 * The dataset layout is @c [x][y][z] in C order. Nonzero values denote
 * solid voxels. The file is opened read-only with the MPI-IO driver.
 * The constructor is collective; every rank throws on failure.
 * Reads use collective MPI-IO transfers, so every rank must call
 * @ref read equally often. A rank without data reads an empty box.
 */
class VoxelMaskFile {
public:
  /**
   * @throw std::runtime_error naming the path or dataset when either is
   * missing or the dataset is not 3D.
   */
  VoxelMaskFile(std::string const &path, std::string const &dataset,
                boost::mpi::communicator const &comm);
  ~VoxelMaskFile();
  VoxelMaskFile(VoxelMaskFile const &) = delete;
  VoxelMaskFile &operator=(VoxelMaskFile const &) = delete;

  /** @brief Extent of the dataset. */
  Utils::Vector3i shape() const;

  /**
   * @brief Collective read of the box [offset, offset + count).
   * A box with a zero count is an empty selection. The data is read in
   * the file's own type, so the transfer needs no type conversion.
   * @return Row-major (x slowest), values converted to 0/1.
   * @throw std::runtime_error if the box lies outside of the dataset,
   * after taking part in the collective transfer with an empty box.
   */
  std::vector<std::uint8_t> read(Utils::Vector3i const &offset,
                                 Utils::Vector3i const &count) const;

  /** @brief Whether the last read used a collective transfer. */
  bool last_read_collective() const;

private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

/** @brief Default read size of @ref read_voxel_mask, in voxels. */
inline constexpr std::size_t voxel_mask_voxels_per_read = std::size_t{1} << 20;

/**
 * @brief Read the stored region of a mask from a file. Collective.
 *
 * Each rank reads its region in boxes of at most @p max_voxels voxels,
 * or of one row along z if a row is longer, and packs them into the bits
 * of @p mask. Ranks with fewer boxes pad with empty reads. Voxels that
 * map to no mask voxel are pore.
 *
 * @throw std::runtime_error on the ranks that failed, after all ranks
 * finished the collective reads.
 */
void read_voxel_mask(VoxelMaskFile const &file, Shapes::VoxelMask &mask,
                     Shapes::VoxelMask::Region const &region,
                     boost::mpi::communicator const &comm,
                     std::size_t max_voxels = voxel_mask_voxels_per_read);

} // namespace IO
