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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE VoxelMaskFile test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "io/voxel_mask/VoxelMaskFile.hpp"

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <highfive/highfive.hpp>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string test_file_path() {
  auto const *tmpdir = std::getenv("TMPDIR");
  auto const dir = std::filesystem::path(tmpdir ? tmpdir : ".");
  return (dir / "VoxelMaskFile_test.h5").string();
}

int expected(int x, int y, int z) { return (x * 100 + y * 10 + z) % 2; }

void write_files(std::string const &path) {
  std::vector<std::size_t> const dims{4, 5, 6};
  std::vector<std::uint8_t> u8;
  std::vector<std::int32_t> i32;
  std::vector<std::int8_t> b8;
  std::vector<double> f64;
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 5; ++y)
      for (int z = 0; z < 6; ++z) {
        auto const v = expected(x, y, z);
        u8.push_back(static_cast<std::uint8_t>(v));
        i32.push_back(v ? 7 : 0);
        b8.push_back(static_cast<std::int8_t>(v));
        // negative zero is pore
        f64.push_back(v ? 0.5 : ((x + y) % 2 ? -0. : 0.));
      }
  HighFive::File file(path, HighFive::File::Truncate);
  file.createDataSet("u8", HighFive::DataSpace(dims),
                     HighFive::create_datatype<std::uint8_t>())
      .write_raw(u8.data());
  file.createDataSet("i32", HighFive::DataSpace(dims),
                     HighFive::create_datatype<std::int32_t>())
      .write_raw(i32.data());
  // h5py-compatible bool: HDF5 enum over int8 with FALSE=0, TRUE=1
  auto const bool_type = HighFive::create_enum_boolean();
  file.createDataSet("flag", HighFive::DataSpace(dims), bool_type)
      .write_raw(b8.data(), bool_type);
  file.createDataSet("f64", HighFive::DataSpace(dims),
                     HighFive::create_datatype<double>())
      .write_raw(f64.data());
}

void check_hyperslab(IO::VoxelMaskFile const &mask) {
  auto const shape = mask.shape();
  BOOST_CHECK_EQUAL(shape[0], 4);
  BOOST_CHECK_EQUAL(shape[1], 5);
  BOOST_CHECK_EQUAL(shape[2], 6);
  auto const values = mask.read({1, 2, 3}, {2, 2, 2});
  BOOST_REQUIRE_EQUAL(values.size(), 8u);
  std::size_t i = 0;
  for (int x = 1; x < 3; ++x)
    for (int y = 2; y < 4; ++y)
      for (int z = 3; z < 5; ++z)
        BOOST_CHECK_EQUAL(static_cast<int>(values[i++]), expected(x, y, z));
}
} // namespace

struct Fixture {
  boost::mpi::communicator comm;
  std::string path = test_file_path();
  Fixture() {
    if (comm.rank() == 0)
      write_files(path);
    comm.barrier();
  }
  ~Fixture() {
    comm.barrier();
    if (comm.rank() == 0)
      std::filesystem::remove(path);
  }
};

BOOST_FIXTURE_TEST_CASE(reads_hyperslab, Fixture) {
  IO::VoxelMaskFile const mask(path, "u8", comm);
  check_hyperslab(mask);
}

BOOST_FIXTURE_TEST_CASE(reads_bool_int32_and_float_datasets, Fixture) {
  check_hyperslab(IO::VoxelMaskFile(path, "i32", comm));
  check_hyperslab(IO::VoxelMaskFile(path, "flag", comm));
  check_hyperslab(IO::VoxelMaskFile(path, "f64", comm));
}

BOOST_FIXTURE_TEST_CASE(unequal_read_counts_are_collective, Fixture) {
  // rank 0 reads three blocks, the other ranks pad with empty selections
  IO::VoxelMaskFile const mask(path, "u8", comm);
  for (int x = 0; x < 3; ++x) {
    if (comm.rank() == 0) {
      auto const values = mask.read({x, 1, 2}, {1, 3, 4});
      BOOST_REQUIRE_EQUAL(values.size(), 12u);
      std::size_t i = 0;
      for (int y = 1; y < 4; ++y)
        for (int z = 2; z < 6; ++z)
          BOOST_CHECK_EQUAL(static_cast<int>(values[i++]), expected(x, y, z));
    } else {
      BOOST_CHECK(mask.read({0, 0, 0}, {0, 0, 0}).empty());
    }
    BOOST_CHECK(mask.last_read_collective());
  }
}

BOOST_FIXTURE_TEST_CASE(invalid_region_throws_after_collective_read, Fixture) {
  IO::VoxelMaskFile const mask(path, "u8", comm);
  if (comm.rank() == 0) {
    BOOST_CHECK_THROW(mask.read({3, 0, 0}, {2, 1, 1}), std::runtime_error);
  } else {
    BOOST_CHECK(mask.read({0, 0, 0}, {1, 1, 1}).size() == 1u);
  }
}

BOOST_FIXTURE_TEST_CASE(missing_dataset_throws, Fixture) {
  BOOST_CHECK_THROW(IO::VoxelMaskFile(path, "nope", comm), std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(missing_file_throws, Fixture) {
  BOOST_CHECK_THROW(IO::VoxelMaskFile(path + ".missing", "u8", comm),
                    std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(reads_regions_with_unequal_read_counts, Fixture) {
  Shapes::VoxelMask::Geometry g;
  g.mask_shape = {4, 5, 6};
  g.voxel_size = 1.;
  g.box_l = {4., 5., 6.};
  g.periodic = {true, true, false};
  g.halo = 1;
  // rank 0: a region that wraps around the periodic x and y faces,
  // rank 1: a region beyond the non-periodic z face, which reads nothing,
  // other ranks: the whole box
  Utils::Vector3d lower{0., 0., 0.};
  Utils::Vector3d upper{4., 5., 6.};
  if (comm.rank() == 0) {
    upper = {2., 3., 6.};
  } else if (comm.rank() == 1) {
    lower = {0., 0., 8.};
    upper = {4., 5., 10.};
  }
  auto const region = Shapes::VoxelMask::stored_region(g, lower, upper);
  IO::VoxelMaskFile const file(path, "u8", comm);
  for (std::size_t max_voxels :
       {std::size_t{1}, std::size_t{7}, std::size_t{40}, std::size_t{1000}}) {
    Shapes::VoxelMask mask;
    mask.set_geometry(g);
    IO::read_voxel_mask(file, mask, region, comm, max_voxels);
    BOOST_REQUIRE_EQUAL(
        mask.local_voxel_count(),
        static_cast<std::size_t>(region.extent[0] * region.extent[1] *
                                 region.extent[2]));
    auto n_solid = 0;
    for (int i = 0; i < region.extent[0]; ++i)
      for (int j = 0; j < region.extent[1]; ++j)
        for (int k = 0; k < region.extent[2]; ++k) {
          Utils::Vector3i const u{region.lower[0] + i, region.lower[1] + j,
                                  region.lower[2] + k};
          Utils::Vector3i m;
          bool pore = false;
          for (int a = 0; a < 3; ++a) {
            m[a] = Shapes::VoxelMask::mask_index(g, a, u[a]);
            pore = pore or m[a] < 0;
          }
          auto const ref = not pore and expected(m[0], m[1], m[2]) == 1;
          Utils::Vector3d const pos{u[0] + 0.5, u[1] + 0.5, u[2] + 0.5};
          BOOST_REQUIRE_EQUAL(mask.is_inside(pos), ref);
          n_solid += ref ? 1 : 0;
        }
    if (comm.rank() == 1) {
      BOOST_CHECK_EQUAL(n_solid, 0);
    } else {
      BOOST_CHECK_GT(n_solid, 0);
    }
  }
}

int main(int argc, char **argv) {
  boost::mpi::environment mpi_env(argc, argv);
  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
