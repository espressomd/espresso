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

#include "io/writer/hdf5_patches.hpp" // must appear first

#include "VoxelMaskFile.hpp"

#include <utils/Vector.hpp>

#include <shapes/VoxelMask.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/operations.hpp>

#include <highfive/highfive.hpp>

#include <H5Dpublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <hdf5.h>
#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace IO {

namespace {
/** Owner of an HDF5 identifier. */
template <herr_t (*Close)(hid_t)> class Handle {
public:
  explicit Handle(hid_t id) : m_id(id) {}
  Handle(Handle const &) = delete;
  Handle &operator=(Handle const &) = delete;
  ~Handle() {
    if (m_id >= 0) {
      Close(m_id);
    }
  }
  hid_t get() const { return m_id; }

private:
  hid_t m_id;
};
} // namespace

struct VoxelMaskFile::Impl {
  HighFive::File file;
  HighFive::DataSet dataset;
  Utils::Vector3i shape;
  H5T_class_t type_class;
  std::size_t type_size;
  /** Native type of the dataset, used as memory type. */
  hid_t mem_type = -1;
  std::string description;
  bool last_read_collective = false;

  Impl(std::string const &path, std::string const &name,
       boost::mpi::communicator const &comm)
      : file(open_file(path, comm)), dataset(open_dataset(path, name)),
        description("dataset '" + name + "' in file '" + path + "'") {
    auto const dims = dataset.getSpace().getDimensions();
    if (dims.size() != 3u) {
      throw std::runtime_error("The " + description + " must be 3D, got " +
                               std::to_string(dims.size()) + " dimensions");
    }
    for (int i = 0; i < 3; ++i) {
      shape[i] = static_cast<int>(dims[i]);
    }
    auto const type = dataset.getDataType();
    type_class = H5Tget_class(type.getId());
    if (type_class != H5T_INTEGER and type_class != H5T_ENUM and
        type_class != H5T_FLOAT) {
      throw std::runtime_error("The " + description +
                               " has an unsupported data type");
    }
    if (type_class == H5T_ENUM and type.getSize() != 1u) {
      throw std::runtime_error("The " + description +
                               " is an enum with more than one byte");
    }
    if (type_class == H5T_FLOAT) {
      // values are compared to zero as float or double
      mem_type = H5Tcopy(type.getSize() == sizeof(float) ? H5T_NATIVE_FLOAT
                                                         : H5T_NATIVE_DOUBLE);
    } else {
      // integers and enums are zero iff all their bytes are zero
      mem_type = H5Tget_native_type(type.getId(), H5T_DIR_ASCEND);
    }
    if (mem_type < 0) {
      throw std::runtime_error("The " + description +
                               " has an unsupported data type");
    }
    type_size = H5Tget_size(mem_type);
  }

  ~Impl() {
    if (mem_type >= 0) {
      H5Tclose(mem_type);
    }
  }
  Impl(Impl const &) = delete;
  Impl &operator=(Impl const &) = delete;

private:
  HighFive::File open_file(std::string const &path,
                           boost::mpi::communicator const &comm) {
    HighFive::FileAccessProps fapl;
    fapl.add(HighFive::MPIOFileAccess(MPI_Comm(comm), MPI_INFO_NULL));
    try {
      // failures are reported via exceptions, not on the HDF5 error stack
      HighFive::SilenceHDF5 const silence{};
      return HighFive::File(path, HighFive::File::ReadOnly, fapl);
    } catch (HighFive::Exception const &err) {
      throw std::runtime_error("Cannot open HDF5 file '" + path +
                               "': " + err.what());
    }
  }

  HighFive::DataSet open_dataset(std::string const &path,
                                 std::string const &name) {
    try {
      HighFive::SilenceHDF5 const silence{};
      return file.getDataSet(name);
    } catch (HighFive::Exception const &err) {
      throw std::runtime_error("Cannot find dataset '" + name + "' in file '" +
                               path + "': " + err.what());
    }
  }
};

VoxelMaskFile::VoxelMaskFile(std::string const &path,
                             std::string const &dataset,
                             boost::mpi::communicator const &comm)
    : m_impl(std::make_unique<Impl>(path, dataset, comm)) {}

VoxelMaskFile::~VoxelMaskFile() = default;

Utils::Vector3i VoxelMaskFile::shape() const { return m_impl->shape; }

bool VoxelMaskFile::last_read_collective() const {
  return m_impl->last_read_collective;
}

template <typename T>
static void floats_to_mask(std::vector<unsigned char> const &raw,
                           std::vector<std::uint8_t> &mask) {
  for (std::size_t i = 0; i < mask.size(); ++i) {
    T value;
    std::memcpy(&value, raw.data() + i * sizeof(T), sizeof(T));
    mask[i] = static_cast<std::uint8_t>(value != T{0});
  }
}

std::vector<std::uint8_t>
VoxelMaskFile::read(Utils::Vector3i const &offset,
                    Utils::Vector3i const &count) const {
  auto &impl = *m_impl;
  std::optional<std::string> error;
  std::size_t total = 1;
  for (int i = 0; i < 3; ++i) {
    if (offset[i] < 0 or count[i] < 0 or offset[i] + count[i] > impl.shape[i]) {
      error = "Requested region lies outside of the " + impl.description;
    }
    total *= static_cast<std::size_t>(std::max(count[i], 0));
  }
  if (error) {
    // take part in the collective transfer with an empty selection
    total = 0;
  }

  HighFive::SilenceHDF5 const silence{};
  Handle<H5Sclose> const file_space(H5Dget_space(impl.dataset.getId()));
  hsize_t const n_elements = total;
  Handle<H5Sclose> const mem_space(H5Screate_simple(1, &n_elements, nullptr));
  Handle<H5Pclose> const xfer(H5Pcreate(H5P_DATASET_XFER));
  auto status = std::min({file_space.get(), mem_space.get(), xfer.get()}) < 0
                    ? herr_t{-1}
                    : H5Pset_dxpl_mpio(xfer.get(), H5FD_MPIO_COLLECTIVE);
  if (status >= 0) {
    if (total == 0u) {
      status = std::min(H5Sselect_none(file_space.get()),
                        H5Sselect_none(mem_space.get()));
    } else {
      hsize_t const start[3] = {static_cast<hsize_t>(offset[0]),
                                static_cast<hsize_t>(offset[1]),
                                static_cast<hsize_t>(offset[2])};
      hsize_t const extent[3] = {static_cast<hsize_t>(count[0]),
                                 static_cast<hsize_t>(count[1]),
                                 static_cast<hsize_t>(count[2])};
      status = H5Sselect_hyperslab(file_space.get(), H5S_SELECT_SET, start,
                                   nullptr, extent, nullptr);
    }
  }
  // raw buffer of the native type; one dummy byte for an empty selection
  std::vector<unsigned char> raw(
      std::max(total * impl.type_size, std::size_t{1}));
  if (status >= 0) {
    // collective: every rank calls H5Dread, a rank without data or with an
    // invalid box selects nothing
    status = H5Dread(impl.dataset.getId(), impl.mem_type, mem_space.get(),
                     file_space.get(), xfer.get(), raw.data());
  }
  H5D_mpio_actual_io_mode_t io_mode = H5D_MPIO_NO_COLLECTIVE;
  if (xfer.get() >= 0) {
    H5Pget_mpio_actual_io_mode(xfer.get(), &io_mode);
  }
  impl.last_read_collective = io_mode != H5D_MPIO_NO_COLLECTIVE;
  if (error) {
    throw std::runtime_error(*error);
  }
  if (status < 0) {
    throw std::runtime_error("Cannot read from the " + impl.description);
  }

  std::vector<std::uint8_t> mask(total);
  if (impl.type_class == H5T_FLOAT and impl.type_size == sizeof(float)) {
    floats_to_mask<float>(raw, mask);
  } else if (impl.type_class == H5T_FLOAT) {
    floats_to_mask<double>(raw, mask);
  } else {
    auto const size = impl.type_size;
    for (std::size_t i = 0; i < total; ++i) {
      auto const *first = raw.data() + i * size;
      mask[i] = static_cast<std::uint8_t>(
          std::any_of(first, first + size, [](auto b) { return b != 0u; }));
    }
  }
  return mask;
}

namespace {
/** Contiguous run of voxels along one axis, mapped to the mask. */
struct Segment {
  int region_start;
  int mask_start;
  int length;
};

/** Split the stored region along @p axis into runs of mask indices. */
std::vector<Segment> segments(Shapes::VoxelMask::Geometry const &g,
                              Shapes::VoxelMask::Region const &region,
                              int axis) {
  std::vector<Segment> result;
  for (int i = 0; i < region.extent[axis]; ++i) {
    auto const m =
        Shapes::VoxelMask::mask_index(g, axis, region.lower[axis] + i);
    if (m < 0) {
      continue;
    }
    if (not result.empty()) {
      auto &last = result.back();
      if (last.region_start + last.length == i and
          last.mask_start + last.length == m) {
        ++last.length;
        continue;
      }
    }
    result.push_back({i, m, 1});
  }
  return result;
}

/** Box to read: position in the mask, extent, position in the region. */
struct Block {
  Utils::Vector3i mask_lower;
  Utils::Vector3i extent;
  Utils::Vector3i region_lower;
};

/** Split the region into boxes of at most @p max_voxels voxels or a row. */
std::vector<Block> blocks(Shapes::VoxelMask::Geometry const &g,
                          Shapes::VoxelMask::Region const &region,
                          std::size_t max_voxels) {
  auto const sx = segments(g, region, 0);
  auto const sy = segments(g, region, 1);
  auto const sz = segments(g, region, 2);
  std::vector<Block> result;
  for (auto const &x : sx) {
    for (auto const &y : sy) {
      for (auto const &z : sz) {
        auto const row = static_cast<std::size_t>(z.length);
        auto const plane = static_cast<std::size_t>(y.length) * row;
        // number of x planes, or of y rows if a plane is too large
        auto const n_x = static_cast<int>(
            std::min(std::max(max_voxels / plane, std::size_t{1}),
                     static_cast<std::size_t>(x.length)));
        auto const n_y =
            (plane <= max_voxels)
                ? y.length
                : static_cast<int>(std::max(max_voxels / row, std::size_t{1}));
        for (int i = 0; i < x.length; i += n_x) {
          for (int j = 0; j < y.length; j += n_y) {
            auto const dx = std::min(n_x, x.length - i);
            auto const dy = std::min(n_y, y.length - j);
            result.push_back(
                {{x.mask_start + i, y.mask_start + j, z.mask_start},
                 {dx, dy, z.length},
                 {x.region_start + i, y.region_start + j, z.region_start}});
          }
        }
      }
    }
  }
  return result;
}
} // namespace

void read_voxel_mask(VoxelMaskFile const &file, Shapes::VoxelMask &mask,
                     Shapes::VoxelMask::Region const &region,
                     boost::mpi::communicator const &comm,
                     std::size_t max_voxels) {
  auto const plan = blocks(mask.geometry(), region, max_voxels);
  auto const n_reads = boost::mpi::all_reduce(
      comm, plan.size(), boost::mpi::maximum<std::size_t>());
  mask.allocate_local_data(region);
  std::exception_ptr error;
  for (std::size_t i = 0; i < n_reads; ++i) {
    if (i >= plan.size() or error) {
      // collective transfer: take part with an empty selection
      try {
        file.read({0, 0, 0}, {0, 0, 0});
      } catch (...) {
        if (not error) {
          error = std::current_exception();
        }
      }
      continue;
    }
    auto const &block = plan[i];
    try {
      auto const voxels = file.read(block.mask_lower, block.extent);
      mask.set_local_block(block.region_lower, block.extent, voxels);
    } catch (...) {
      error = std::current_exception();
    }
  }
  if (error) {
    std::rethrow_exception(error);
  }
}

} // namespace IO
