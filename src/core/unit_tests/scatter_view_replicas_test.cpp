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

#define BOOST_TEST_MODULE "ScatterView replica reduction"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cell_system/CellStructure.hpp"
#include "cell_system/scatter_view_replicas.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_ScatterView.hpp>

#include <cstddef>
#include <type_traits>

struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

using ForceType = CellStructure::ForceType;
using ScatterForce = CellStructure::ScatterForce;
using execution_space = CellStructure::execution_space;

/**
 * True if the force scatter view must take the fused reduction. With the
 * OpenMP backend, Kokkos duplicates the scatter view; if a Kokkos update
 * changes the ScatterView template signature or the replica storage, the
 * fused path silently falls back to contribute() + reset() and this test
 * fails instead.
 */
static constexpr bool fused_path_expected() {
#ifdef KOKKOS_ENABLE_OPENMP
  return std::is_same_v<execution_space, Kokkos::OpenMP>;
#else
  return false;
#endif
}

BOOST_AUTO_TEST_CASE(force_scatter_view_type) {
  if constexpr (fused_path_expected()) {
    static_assert(ScatterViewReplicaReduction::is_duplicated_sum_scatter_view<
                  ScatterForce>::value);
  }
}

/** Scatter integer-valued contributions from all threads. */
static void scatter(ScatterForce &scatter_view, std::size_t n) {
  Kokkos::parallel_for("scatter",
                       Kokkos::RangePolicy<execution_space>(0, 4 * n),
                       [=](std::size_t const k) {
                         auto access = scatter_view.access();
                         auto const i = k % n;
                         for (int d = 0; d < 3; ++d) {
                           access(i, d) += static_cast<double>(k % 7 + d);
                         }
                       });
  execution_space{}.fence();
}

BOOST_AUTO_TEST_CASE(fused_reduction) {
  // block of the reduction is 512 values, i.e. 170.67 particles
  for (std::size_t const n : {1ul, 170ul, 171ul, 512ul, 1000ul, 4099ul}) {
    ForceType dest("dest", n);
    ForceType reference("reference", n);
    Kokkos::deep_copy(dest, 1.);
    Kokkos::deep_copy(reference, 1.);
    ScatterForce scatter_dest(dest);
    ScatterForce scatter_reference(reference);
    scatter(scatter_dest, n);
    scatter(scatter_reference, n);

    auto const fused =
        ScatterViewReplicaReduction::contribute_and_reset_replicas(scatter_dest,
                                                                   dest);
    Kokkos::Experimental::contribute(reference, scatter_reference);

    BOOST_CHECK_EQUAL(fused, fused_path_expected());
    if (fused) {
      BOOST_CHECK(ScatterViewReplicaReduction::replicas_are_zero(scatter_dest));
    }
    for (std::size_t i = 0u; i < n; ++i) {
      for (int d = 0; d < 3; ++d) {
        BOOST_REQUIRE_EQUAL(dest(i, d), reference(i, d));
      }
    }

    // a second round on the zeroed replicas, without reset()
    if (fused) {
      scatter(scatter_dest, n);
      scatter_reference.reset();
      scatter(scatter_reference, n);
      ScatterViewReplicaReduction::contribute_and_reset_replicas(scatter_dest,
                                                                 dest);
      Kokkos::Experimental::contribute(reference, scatter_reference);
      for (std::size_t i = 0u; i < n; ++i) {
        for (int d = 0; d < 3; ++d) {
          BOOST_REQUIRE_EQUAL(dest(i, d), reference(i, d));
        }
      }
    }
  }
}
