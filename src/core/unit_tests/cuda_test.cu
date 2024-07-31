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

#define BOOST_TEST_MODULE "CUDA interface tests"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "cuda/init.hpp"
#include "cuda/utils.cuh"
#include "cuda/utils.hpp"
#include "errorhandling.hpp"

#include "cuda/CudaHostAllocator.hpp"

#include <cuda.h>

#include <boost/test/unit_test.hpp>

#include <cstddef>
#include <initializer_list>
#include <optional>
#include <string>

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

std::optional<std::string> read_pending_cuda_errors() {
  auto const CU_err = cudaGetLastError();
  if (CU_err != cudaSuccess) {
    auto const message = std::string(cudaGetErrorString(CU_err));
    return {"There is a pending CUDA error: \"" + message + "\""};
  }
  return std::nullopt;
}

void setup() {}
void teardown() {
  auto error = read_pending_cuda_errors();
  BOOST_REQUIRE_MESSAGE(not error.has_value(), error.value_or(""));
}

namespace Testing::non_sticky_cuda_error {

/** @brief Trigger a non-sticky CUDA error for testing purposes. */
void trigger() { cudaSetDevice(-1); }

/** @brief Clear a non-sticky CUDA error raised by @ref trigger. */
void clear() {
  auto const error_code = cudaGetLastError();
  BOOST_REQUIRE_MESSAGE(error_code == cudaSuccess or
                            error_code == cudaErrorInvalidDevice,
                        "An unexpected CUDA error was pending!");
}

} // namespace Testing::non_sticky_cuda_error

#ifdef P3M
dim3 p3m_make_grid(unsigned int n_blocks);
#endif

static auto fixture = boost::unit_test::fixture(&setup, &teardown);

BOOST_AUTO_TEST_SUITE(suite, *boost::unit_test::precondition(has_gpu))

BOOST_AUTO_TEST_CASE(gpu_fixture, *fixture) {
  auto error1 = read_pending_cuda_errors();
  BOOST_REQUIRE(not error1.has_value());

  // check we can raise and clear non-sticky CUDA errors
  Testing::non_sticky_cuda_error::trigger();
  Testing::non_sticky_cuda_error::clear();
  auto error2 = read_pending_cuda_errors();
  BOOST_REQUIRE(not error2.has_value());

  // check fixture can handle the default non-sticky CUDA error
  Testing::non_sticky_cuda_error::trigger();
  auto ref_what3{"There is a pending CUDA error: \"invalid device ordinal\""};
  auto error3 = read_pending_cuda_errors();
  BOOST_REQUIRE(error3.has_value());
  BOOST_REQUIRE_EQUAL(error3.value(), ref_what3);
  // sticky error should have been cleared
  error3 = read_pending_cuda_errors();
  BOOST_REQUIRE(not error3.has_value());

  // check fixture can handle a custom non-sticky CUDA error
  cudaMallocHost(nullptr, std::size_t(0u));
  auto ref_what4{"There is a pending CUDA error: \"invalid argument\""};
  auto error4 = read_pending_cuda_errors();
  BOOST_REQUIRE(error4.has_value());
  BOOST_REQUIRE_EQUAL(error4.value(), ref_what4);
  // sticky error should have been cleared
  error4 = read_pending_cuda_errors();
  BOOST_REQUIRE(not error4.has_value());
}

static int fatal_error_counter = 0;
static void increment_counter() noexcept { ++fatal_error_counter; }

BOOST_AUTO_TEST_CASE(gpu_interface, *fixture) {
  fatal_error_counter = 0;
  auto local_error_counter = 0;

  try {
    throw cuda_fatal_error("message");
  } catch (cuda_fatal_error &err) {
    std::string const what = "message";
    BOOST_CHECK_EQUAL(err.what(), what);
    BOOST_CHECK_EQUAL(err.get_terminate(), &errexit);
    err.set_terminate(nullptr);
    BOOST_CHECK_EQUAL(err.get_terminate(), nullptr);
    err.set_terminate(increment_counter);
    BOOST_CHECK_EQUAL(err.get_terminate(), &increment_counter);
    BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
  }
  ++local_error_counter;
  BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);

  // -----------------------

  auto error_caught = false;
  auto const block = dim3{1, 2, 3};
  auto const grid = dim3{4, 5, 6};
  cuda_check_errors_exit(block, grid, "", "", 0u); // should not throw
  try {
    Testing::non_sticky_cuda_error::trigger();
    // should clear the CUDA error flag and throw a fatal error
    cuda_check_errors_exit(block, grid, "cudaSetDevice()", "filename.cu", 4u);
  } catch (cuda_fatal_error &err) {
    error_caught = true;
    err.set_terminate(increment_counter);
    std::string const what =
        "CUDA error: \"invalid device ordinal\" while calling "
        "cudaSetDevice() with block: <1,2,3>, grid: <4,5,6> in filename.cu:4";
    BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
    BOOST_CHECK_EQUAL(err.what(), what);
    BOOST_CHECK_EQUAL(cudaGetLastError(), cudaSuccess);
  }
  ++local_error_counter;
  BOOST_REQUIRE(error_caught);
  BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);

  // -----------------------

  error_caught = false;
  cuda_safe_mem_exit(cudaSuccess, "", 0u); // should not throw
  try {
    Testing::non_sticky_cuda_error::trigger();
    cuda_safe_mem_exit(cudaSuccess, "filename.cu", 4u); // should throw
  } catch (cuda_fatal_error &err) {
    error_caught = true;
    err.set_terminate(increment_counter);
    std::string const what =
        "CUDA error: \"invalid device ordinal\" in filename.cu:4. Error "
        "found during memory operation. Possibly however from a failed "
        "operation before the memory operation";
    BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
    BOOST_CHECK_EQUAL(err.what(), what);
  }
  ++local_error_counter;
  BOOST_REQUIRE(error_caught);
  BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);

  // -----------------------

  error_caught = false;
  try {
    cuda_safe_mem_exit(cudaErrorNotPermitted, "filename.cu", 4u);
  } catch (cuda_fatal_error &err) {
    error_caught = true;
    err.set_terminate(increment_counter);
    std::string const what = "CUDA error: \"operation not permitted\" during "
                             "memory operation in filename.cu:4";
    BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
    BOOST_CHECK_EQUAL(err.what(), what);
  }
  ++local_error_counter;
  BOOST_REQUIRE(error_caught);
  BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);

  // -----------------------

  error_caught = false;
  try {
    cuda_safe_mem_exit(cudaErrorInvalidValue, "function_name()", 4u);
  } catch (cuda_fatal_error &err) {
    error_caught = true;
    err.set_terminate(increment_counter);
    std::string const what =
        "CUDA error: \"invalid argument\" during memory operation in "
        "function_name():4. You may have tried to allocate zero memory";
    BOOST_CHECK_EQUAL(fatal_error_counter, local_error_counter);
    BOOST_CHECK_EQUAL(err.what(), what);
  }
  ++local_error_counter;
  BOOST_REQUIRE(error_caught);
  BOOST_REQUIRE_EQUAL(fatal_error_counter, local_error_counter);

  // -----------------------

  error_caught = false;
  BOOST_REQUIRE_EQUAL(stream[0], nullptr);
  cuda_init(); // allocate
  BOOST_REQUIRE_NE(stream[0], nullptr);
  cuda_set_device(0); // reallocate, may or may not result in the same pointer
  BOOST_REQUIRE_NE(stream[0], nullptr);
  auto const old_stream = stream[0];
  try {
    cuda_set_device(-1); // fail to reallocate, pointer remains the same
  } catch (cuda_runtime_error_cuda const &err) {
    error_caught = true;
    std::string const what = "CUDA error: invalid device ordinal";
    BOOST_CHECK_EQUAL(err.what(), what);
  }
  BOOST_REQUIRE(error_caught);
  BOOST_REQUIRE_EQUAL(stream[0], old_stream);

  // -----------------------

  BOOST_REQUIRE_GE(cuda_get_n_gpus(), 1);
  char gpu_name_buffer[260] = {'\0'};
  cuda_get_gpu_name(0, gpu_name_buffer);
  for (int i = 255; i < 260; ++i) {
    BOOST_REQUIRE_EQUAL(gpu_name_buffer[i], '\0');
  }
}

#ifdef P3M

BOOST_AUTO_TEST_CASE(p3m_reshape_grid_test, *fixture) {
  auto constexpr optimal_size = 65536u;

  for (auto cao = 1u; cao <= 3u; ++cao) {
    auto const n_blocks = cao * optimal_size;
    auto const grid = p3m_make_grid(n_blocks);
    BOOST_CHECK_EQUAL(grid.x, optimal_size);
    BOOST_CHECK_EQUAL(grid.y, cao);
    BOOST_CHECK_EQUAL(grid.z, 1u);
  }

  for (auto mul : {2u, 3u, 6u, 12u}) {
    auto const n_blocks = mul * optimal_size + 1u;
    auto const grid = p3m_make_grid(n_blocks);
    BOOST_CHECK_EQUAL(grid.x, n_blocks / (mul + 1u) + ((mul == 2u) ? 0u : 1u));
    BOOST_CHECK_EQUAL(grid.y, (mul + 1u));
    BOOST_CHECK_EQUAL(grid.z, 1u);
  }
}

#endif

BOOST_AUTO_TEST_SUITE_END()
