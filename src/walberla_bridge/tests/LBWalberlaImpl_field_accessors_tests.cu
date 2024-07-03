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
#define BOOST_TEST_MODULE "field accessors test"
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_NO_MAIN

#if defined(__NVCC__)
// Fix for https://i10git.cs.fau.de/walberla/walberla/-/issues/244
#pragma nv_diagnostic push
#pragma nv_diag_suppress 554 // unreachable conversion operator
#endif

#include <boost/test/data/monomorphic.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "../src/lattice_boltzmann/LBWalberlaImpl.hpp"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsDoublePrecision.h"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsDoublePrecisionCUDA.cuh"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsSinglePrecision.h"
#include "../src/lattice_boltzmann/generated_kernels/FieldAccessorsSinglePrecisionCUDA.cuh"

#include <walberla_bridge/Architecture.hpp>
#include <walberla_bridge/BlockAndCell.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpl/list.hpp>

#include <mpi.h>

#include <cuda.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <numeric>
#include <ranges>
#include <span>
#include <sstream>
#include <vector>

static Utils::Vector3i mpi_shape;

boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id) {
  bool has_compatible_device = false;
  int n_devices = 0;
  cudaGetDeviceCount(&n_devices);
  if (n_devices > 0) {
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    if (prop.major >= 3) {
      has_compatible_device = true;
    }
  }
  return has_compatible_device;
}

template <std::ranges::contiguous_range R>
  requires std::ranges::sized_range<R>
boost::test_tools::predicate_result almost_equal(R const &val, R const &ref,
                                                 typename R::value_type atol) {
  assert(val.size() == ref.size());
  assert(not val.empty());
  auto const print_first_n = [](R const &v) {
    auto constexpr n = std::size_t{6u};
    std::stringstream stream;
    stream << "[";
    for (auto i = 0ul, end = std::min(v.size(), n); i < end; ++i) {
      if (i) {
        stream << ", ";
      }
      stream << v[i];
    }
    if (v.size() > n) {
      stream << ", ...";
    }
    stream << "]";
    return stream.str();
  };
  boost::test_tools::predicate_result res(true);
  for (auto i = 0ul; i < val.size(); ++i) {
    if (auto const diff = std::abs(val[i] - ref[i]); diff > atol) {
      res = false;
      res.message() << "val{" << print_first_n(val) << "} and " << "ref{"
                    << print_first_n(ref) << "} mismatch: " << "val[" << i
                    << "]{" << val[i] << "} != " << "ref[" << i << "]{"
                    << ref[i] << "} " << "(difference{" << diff << "} > delta{"
                    << atol << "})";
      break;
    }
  }
  return res;
}

template <typename T>
boost::test_tools::predicate_result almost_equal(walberla::Vector3<T> const &a,
                                                 walberla::Vector3<T> const &b,
                                                 T atol) {
  return almost_equal(std::span(a.data(), 3ul), std::span(b.data(), 3ul), atol);
}

template <typename T>
boost::test_tools::predicate_result almost_equal(walberla::Matrix3<T> const &a,
                                                 walberla::Matrix3<T> const &b,
                                                 T atol) {
  return almost_equal(std::span(a.data(), 9ul), std::span(b.data(), 9ul), atol);
}

template <typename T>
  requires std::floating_point<T>
boost::test_tools::predicate_result almost_equal(T a, T b, T atol) {
  return almost_equal(std::span(&a, 1ul), std::span(&b, 1ul), atol);
}

template <typename FT, lbmpy::Arch Architecture>
class LBWalberlaImplTest : public walberla::LBWalberlaImpl<FT, Architecture> {
public:
  using Base = walberla::LBWalberlaImpl<FT, Architecture>;
  using Base::Base;
  using Base::m_density;
  using Base::m_last_applied_force_field_id;
  using Base::m_pdf_field_id;
  using Base::m_velocity_field_id;
};

template <typename FT, lbmpy::Arch Architecture> struct Fixture {
  std::shared_ptr<::LatticeWalberla> lattice;
  std::shared_ptr<LBWalberlaImplTest<FT, Architecture>> lbfluid;

  Fixture() {
    auto const grid_dim = Utils::Vector3i::broadcast(4);
    auto const viscosity = FT(1.5);
    auto const density = FT(0.9);
    lattice = std::make_shared<::LatticeWalberla>(grid_dim, mpi_shape, 1u);
    lbfluid = std::make_shared<LBWalberlaImplTest<FT, Architecture>>(
        lattice, viscosity, density);
  }

  void runTest() {
    auto const bc = walberla::get_block_and_cell(*lattice, {0, 0, 0}, false);
    BOOST_CHECK(bc);
    auto const cell = bc->cell;
    auto const ci = walberla::CellInterval(cell, cell);
    auto &block = *(bc->block);
    check_getters_setters(block, cell);
    check_getters_setters(block, ci);
  }

  template <typename CellIt>
  void check_getters_setters(walberla::IBlock &block, CellIt const &it) {
    using namespace walberla;
    using PdfField = LBWalberlaImplTest<FT, Architecture>::PdfField;
    using VectorField = LBWalberlaImplTest<FT, Architecture>::VectorField;

    auto constexpr is_interval = std::is_same_v<CellIt, CellInterval>;
    auto constexpr epsilon = std::numeric_limits<FT>::epsilon();
    auto constexpr exact = FT{0};

    auto const make_ref_vector = [](std::initializer_list<FT> values) {
      if constexpr (is_interval) {
        assert(values.size() % 3ul == 0ul);
        return std::vector<FT>(values);
      } else {
        assert(values.size() == 3ul);
        return Vector3<FT>(std::data(values));
      }
    };
    auto const make_ref_matrix = [](std::initializer_list<FT> values) {
      if constexpr (is_interval) {
        assert(values.size() % 9ul == 0ul);
        return std::vector<FT>(values);
      } else {
        assert(values.size() == 9ul);
        return Matrix3<FT>(std::data(values));
      }
    };
    auto const to_number = [](auto const &value) {
      if constexpr (std::is_same_v<decltype(value), FT const &>) {
        return value;
      } else {
        assert(value.size() == 1ul);
        return value[0u];
      }
    };

    auto const density = lbfluid->m_density;
    auto pdf_field = block.template getData<PdfField>(lbfluid->m_pdf_field_id);
    auto velocity_field =
        block.template getData<VectorField>(lbfluid->m_velocity_field_id);
    auto force_field = block.template getData<VectorField>(
        lbfluid->m_last_applied_force_field_id);

    std::conditional_t<is_interval, std::vector<FT>, std::array<FT, 19u>>
        cur_pop = lbm::accessor::Population::get(pdf_field, it);
    std::conditional_t<is_interval, std::vector<FT>, Vector3<FT>> cur_vel;
    std::conditional_t<is_interval, std::vector<FT>, Vector3<FT>> cur_laf;

    auto const reset_fields = [&, initial_pop = cur_pop]() {
      auto const null_vector = Vector3<FT>(FT{0});
      std::array<FT, 19u> pop;
      for (std::size_t i = 0ul; i < pop.size(); ++i) {
        pop[i] = initial_pop[i];
      }
      lbm::accessor::Population::initialize(pdf_field, pop);
      lbm::accessor::Vector::initialize(force_field, null_vector);
      lbm::accessor::Vector::initialize(velocity_field, null_vector);
    };

    {
      auto diag = FT{0};
      auto const zero = FT{0};
      auto const old_pop = lbm::accessor::Population::get(pdf_field, it);
      auto const old_pre = lbm::accessor::PressureTensor::get(pdf_field, it);
      auto const old_laf = lbm::accessor::Vector::get(force_field, it);
      auto const old_rho = lbm::accessor::Density::get(pdf_field, it);
      auto ref_pop = old_pop;
      std::transform(old_pop.begin(), old_pop.end(), ref_pop.begin(),
                     [](auto const &f) { return FT{2} * f; });
      lbm::accessor::Population::set(pdf_field, ref_pop, it);
      auto const new_pop = lbm::accessor::Population::get(pdf_field, it);
      auto const new_pre = lbm::accessor::PressureTensor::get(pdf_field, it);
      auto const new_rho = lbm::accessor::Density::get(pdf_field, it);
      BOOST_CHECK(almost_equal(new_pop, ref_pop, exact));
      // clang-format off
      diag = density * (FT{1} / FT{3});
      auto const old_pre_ref = make_ref_matrix({diag, zero, zero,
                                                zero, diag, zero,
                                                zero, zero, diag});
      BOOST_CHECK(almost_equal(old_pre, old_pre_ref, epsilon));
      diag = density * (FT{2} / FT{3});
      auto const new_pre_ref = make_ref_matrix({diag, zero, zero,
                                                zero, diag, zero,
                                                zero, zero, diag});
      BOOST_CHECK(almost_equal(new_pre, new_pre_ref, epsilon));
      // clang-format on
      auto const old_laf_ref = make_ref_vector({FT{0}, FT{0}, FT{0}});
      auto const new_laf_ref = make_ref_vector({FT{2}, FT{3}, FT{4}});
      lbm::accessor::Vector::set(force_field, new_laf_ref, it);
      auto const new_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(old_laf, old_laf_ref, exact));
      BOOST_CHECK(almost_equal(new_laf, new_laf_ref, exact));
      lbm::accessor::Vector::set(force_field, old_laf_ref, it);
      cur_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(cur_laf, old_laf_ref, exact));
      lbm::accessor::Vector::add_to_all(force_field,
                                        Vector3<FT>(new_laf_ref.data()));
      cur_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(cur_laf, new_laf_ref, exact));
      lbm::accessor::Vector::set(force_field, old_laf_ref, it);
      cur_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(cur_laf, old_laf_ref, exact));
      lbm::accessor::Density::set(pdf_field, {FT{7} * density}, it);
      auto const cur_rho = lbm::accessor::Density::get(pdf_field, it);
      BOOST_CHECK(
          almost_equal(to_number(old_rho), density * FT{1}, FT{20} * epsilon));
      BOOST_CHECK(
          almost_equal(to_number(new_rho), density * FT{2}, FT{20} * epsilon));
      BOOST_CHECK(
          almost_equal(to_number(cur_rho), density * FT{7}, FT{20} * epsilon));
    }
    reset_fields();
    {
      // update cached velocities and recalculate populations in a single step
      auto const old_pop = lbm::accessor::Population::get(pdf_field, it);
      auto const old_vel = lbm::accessor::Vector::get(velocity_field, it);
      auto const ref_vel = make_ref_vector({FT(0.001), FT(0.002), FT(0.003)});
      lbm::accessor::Velocity::set(pdf_field, velocity_field, force_field,
                                   ref_vel, it);
      auto const new_pop = lbm::accessor::Population::get(pdf_field, it);
      auto const new_vel = lbm::accessor::Vector::get(velocity_field, it);
      BOOST_CHECK(almost_equal(new_vel, ref_vel, epsilon));
      auto const new_mom =
          lbm::accessor::MomentumDensity::reduce(pdf_field, force_field);
      auto const ref_mom = Vector3<FT>(
          ref_vel[0u] * density, ref_vel[1u] * density, ref_vel[2u] * density);
      BOOST_CHECK(almost_equal(new_mom, ref_mom, FT{20} * epsilon));
      // update populations and recalculate cached velocities in a single step
      lbm::accessor::Population::set(pdf_field, velocity_field, force_field,
                                     old_pop, it);
      cur_pop = lbm::accessor::Population::get(pdf_field, it);
      cur_vel = lbm::accessor::Vector::get(velocity_field, it);
      BOOST_CHECK(almost_equal(cur_pop, old_pop, exact));
      BOOST_CHECK(almost_equal(cur_vel, old_vel, epsilon));
      cur_vel = lbm::accessor::Velocity::get(pdf_field, force_field, it);
      BOOST_CHECK(almost_equal(cur_vel, old_vel, epsilon));
      lbm::accessor::Population::set(pdf_field, velocity_field, force_field,
                                     new_pop, it);
      cur_pop = lbm::accessor::Population::get(pdf_field, it);
      cur_vel = lbm::accessor::Vector::get(velocity_field, it);
      BOOST_CHECK(almost_equal(cur_pop, new_pop, exact));
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
      cur_vel = lbm::accessor::Velocity::get(pdf_field, force_field, it);
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
      // update forces and recalculate cached velocities in a single step
      auto const ref_laf = make_ref_vector({ref_vel[0u] * FT{2} * density,
                                            ref_vel[1u] * FT{2} * density,
                                            ref_vel[2u] * FT{2} * density});
      lbm::accessor::Population::set(pdf_field, velocity_field, force_field,
                                     old_pop, it);
      lbm::accessor::Force::set(pdf_field, velocity_field, force_field, ref_laf,
                                it);
      cur_pop = lbm::accessor::Population::get(pdf_field, it);
      cur_vel = lbm::accessor::Vector::get(velocity_field, it);
      cur_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(cur_pop, old_pop, exact));
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
      BOOST_CHECK(almost_equal(cur_laf, ref_laf, epsilon));
      cur_vel = lbm::accessor::Velocity::get(pdf_field, force_field, it);
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
      // update velocities and recalculate populations in a single step
      lbm::accessor::Population::set(pdf_field, velocity_field, force_field,
                                     old_pop, it);
      lbm::accessor::Velocity::set(pdf_field, velocity_field, force_field,
                                   ref_vel, it);
      cur_pop = lbm::accessor::Population::get(pdf_field, it);
      cur_vel = lbm::accessor::Vector::get(velocity_field, it);
      cur_laf = lbm::accessor::Vector::get(force_field, it);
      BOOST_CHECK(almost_equal(cur_pop, old_pop, epsilon));
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
      BOOST_CHECK(almost_equal(cur_laf, ref_laf, epsilon));
      cur_vel = lbm::accessor::Velocity::get(pdf_field, force_field, it);
      BOOST_CHECK(almost_equal(cur_vel, new_vel, epsilon));
    }
    reset_fields();
  }
};

BOOST_AUTO_TEST_SUITE(suite, *boost::unit_test::precondition(has_gpu))

BOOST_AUTO_TEST_CASE(test_custom_predicate) {
  std::vector<int> const val = {0, 1, 2, 3, 4, 5, 99, 2};
  std::vector<int> const ref = {0, 1, 2, 3, 4, 5, 7, 80};
  auto const is_true = almost_equal(ref, ref, 0);
  auto const is_false = almost_equal(val, ref, 1);
  BOOST_REQUIRE(is_true);
  BOOST_REQUIRE(not is_false);
  BOOST_CHECK_EQUAL(is_true.message(), "");
  BOOST_REQUIRE_EQUAL(
      is_false.message(),
      "val{[0, 1, 2, 3, 4, 5, ...]} and ref{[0, 1, 2, 3, 4, 5, ...]} "
      "mismatch: val[6]{99} != ref[6]{7} (difference{92} > delta{1})");
}

using test_types = boost::mpl::list<float, double>;

BOOST_AUTO_TEST_CASE_TEMPLATE(macroscopic_accessors, FT, test_types) {
  Fixture<FT, lbmpy::Arch::CPU>().runTest();
  Fixture<FT, lbmpy::Arch::GPU>().runTest();
}

BOOST_AUTO_TEST_SUITE_END()

int main(int argc, char **argv) {
  int n_nodes;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &n_nodes);
  MPI_Dims_create(n_nodes, 3, mpi_shape.data());
  walberla::mpi_init();

  auto const res = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
  MPI_Finalize();
  return res;
}

#if defined(__NVCC__)
#pragma nv_diagnostic pop
#endif
