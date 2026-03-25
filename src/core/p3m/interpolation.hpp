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

#include <config/config.hpp>

#include <utils/index.hpp>
#include <utils/math/bspline.hpp>

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
#include <Kokkos_Core.hpp>
#endif

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

/**
 * @brief Interpolation weights for one point.
 *
 * Interpolation weights and grid offset for one point.
 *
 * @tparam cao Interpolation order.
 */
template <int cao> struct InterpolationWeights {
  /** Linear index of the corner of the interpolation cube. */
  int ind;
  /** Weights for the directions */
  Utils::Array<double, cao> w_x, w_y, w_z;
};

template <int cao> struct InterpolationWeightsView {
  /** Linear index of the corner of the interpolation cube. */
  int ind;
  /** Weights for the directions */
  std::span<double const, cao> w_x, w_y, w_z;
};

/**
 * @brief Cache for interpolation weights.
 *
 * This is a storage container for interpolation weights of
 * type InterpolationWeights.
 */
class p3m_interpolation_cache {
  int m_cao = 0;
  /** Charge fractions for mesh assignment. */
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  Kokkos::View<double *, Kokkos::LayoutRight, Kokkos::HostSpace> ca_frac;
  std::vector<double> ca_frac_spillover;
#else
  std::vector<double> ca_frac;
#endif
  /** index of first mesh point for charge assignment. */
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  Kokkos::View<int *, Kokkos::LayoutRight, Kokkos::HostSpace> ca_fmp;
  std::vector<int> ca_fmp_spillover;
#else
  std::vector<int> ca_fmp;
#endif

public:
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
  p3m_interpolation_cache()
      : ca_frac("p3m_interpolation_cache::ca_frac", 0, 0),
        ca_fmp("p3m_interpolation_cache::ca_fmp", 0, 0) {}
#endif

  /**
   * @brief Number of points in the cache.
   * @return Number of points currently in the cache.
   */
  auto size() const {
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    return ca_fmp.size() + ca_fmp_spillover.size();
#else
    return ca_fmp.size();
#endif
  }

  /**
   * @brief Charge assignment order the weights are for.
   * @return The charge assignment order.
   */
  auto cao() const { return m_cao; }

  /**
   * @brief Fill cache with zero-initialized data.
   * Meant for the parallel weight interpolation calculation.
   * @param size Number of elements.
   */
  void zfill(std::size_t size) {
    auto const size_frac = size * static_cast<std::size_t>(m_cao * 3);
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    if (ca_fmp.size() != size or ca_frac.size() != size_frac) {
      Kokkos::realloc(Kokkos::WithoutInitializing, ca_fmp, size);
      Kokkos::realloc(Kokkos::WithoutInitializing, ca_frac, size_frac);
    }
    assert(ca_frac_spillover.empty());
    assert(ca_fmp_spillover.empty());
    // data must be contiguous in memory
    assert(size == 0 or
           (std::distance(&ca_fmp(0), &ca_fmp(size - 1)) == size - 1));
#else
    assert(ca_frac.empty());
    assert(ca_fmp.empty());
    ca_fmp.resize(size);
    ca_frac.resize(size_frac);
#endif
  }

  /**
   * @brief Push back weights for one point.
   *
   * @tparam cao Interpolation order has to match the order
   *         set at last call to @ref p3m_interpolation_cache::reset.
   * @param weights Interpolation weights to store.
   */
  template <int cao> void store(InterpolationWeights<cao> const &weights) {
    assert(cao == m_cao);

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    ca_fmp_spillover.emplace_back(weights.ind);
    auto it = std::back_inserter(ca_frac_spillover);
    std::ranges::copy(weights.w_x, it);
    std::ranges::copy(weights.w_y, it);
    std::ranges::copy(weights.w_z, it);
#else
    ca_fmp.emplace_back(weights.ind);
    auto it = std::back_inserter(ca_frac);
    std::ranges::copy(weights.w_x, it);
    std::ranges::copy(weights.w_y, it);
    std::ranges::copy(weights.w_z, it);
#endif
  }

  /**
   * @brief Insert weights for one point.
   *
   * @tparam cao Interpolation order has to match the order
   *         set at last call to @ref p3m_interpolation_cache::reset.
   * @param p_index Particle index.
   * @param weights Interpolation weights to store.
   */
  template <int cao>
  void store_at(std::size_t p_index, InterpolationWeights<cao> const &weights) {
    assert(cao == m_cao);
    assert(p_index < ca_fmp.size());

    auto copy_weights = [](InterpolationWeights<cao> const &src, auto dst) {
      std::copy_n(src.w_x.data(), cao, dst);
      std::advance(dst, cao);
      std::copy_n(src.w_y.data(), cao, dst);
      std::advance(dst, cao);
      std::copy_n(src.w_z.data(), cao, dst);
    };

#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    ca_fmp(p_index) = weights.ind;
#else
    ca_fmp[p_index] = weights.ind;
#endif
    auto const offset = p_index * static_cast<std::size_t>(cao * 3);
    copy_weights(weights, ca_frac.data() + offset);
  }

  /**
   * @brief Load entry from the cache.
   *
   * This loads an entry at an index from the cache,
   * the entries are indexed by the order they were stored.
   *
   * @tparam cao Interpolation order has to match the order
   *         set at last call to @ref p3m_interpolation_cache::reset.
   * @param p_index Index of the entry to load.
   * @return Interpolation weights.
   */
  template <int cao>
  InterpolationWeightsView<cao> load(std::size_t p_index) const {
    assert(cao == m_cao);
    assert(p_index < size());

    int index = -1;
    double const *ptr = nullptr;
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    if (p_index >= ca_fmp.size()) {
      p_index -= ca_fmp.size();
      index = ca_fmp_spillover[p_index];
      auto const offset = p_index * static_cast<std::size_t>(cao * 3);
      ptr = std::as_const(ca_frac_spillover).data() + offset;
    } else {
      index = ca_fmp(p_index);
      auto const offset = p_index * static_cast<std::size_t>(cao * 3);
      ptr = std::as_const(ca_frac).data() + offset;
    }
#else
    index = ca_fmp[p_index];
    auto const offset = p_index * static_cast<std::size_t>(cao * 3);
    ptr = std::as_const(ca_frac).data() + offset;
#endif
    std::span<double const, cao> w_x(ptr, cao);
    std::advance(ptr, cao);
    std::span<double const, cao> w_y(ptr, cao);
    std::advance(ptr, cao);
    std::span<double const, cao> w_z(ptr, cao);
    return {index, w_x, w_y, w_z};
  }

  /**
   * @brief Reset the cache.
   *
   * @param cao Interpolation order.
   */
  void reset(int cao) {
    m_cao = cao;
#ifdef ESPRESSO_SHARED_MEMORY_PARALLELISM
    Kokkos::realloc(Kokkos::WithoutInitializing, ca_fmp, 0ul);
    Kokkos::realloc(Kokkos::WithoutInitializing, ca_frac, 0ul);
    ca_frac_spillover.clear();
    ca_fmp_spillover.clear();
#else
    ca_frac.clear();
    ca_fmp.clear();
#endif
  }
};

/**
 * @brief Calculate the P-th order interpolation weights.
 *
 * As described in from @cite hockney88a 5-189 (or 8-61).
 * The weights are also tabulated in @cite deserno98a @cite deserno98b.
 */
template <int cao, Utils::MemoryOrder memory_order>
InterpolationWeights<cao>
p3m_calculate_interpolation_weights(const Utils::Vector3d &position,
                                    const Utils::Vector3d &ai,
                                    P3MLocalMesh const &local_mesh) {
  /** position shift for calc. of first assignment mesh point. */
  auto const pos_shift = std::floor((cao - 1) / 2.0) - (cao % 2) / 2.0;

  /* distance to nearest mesh point */
  Utils::Vector3d dist;

  /* nearest mesh point */
  Utils::Vector3i nmp;

  for (unsigned int d = 0; d < 3; d++) {
    /* particle position in mesh coordinates */
    auto const pos = ((position[d] - local_mesh.ld_pos[d]) * ai[d]) - pos_shift;

    nmp[d] = static_cast<int>(pos);

    /* distance to nearest mesh point */
    dist[d] = (pos - nmp[d]) - 0.5;
  }

  InterpolationWeights<cao> ret;

  /* 3d-array index of nearest mesh point */
  ret.ind = Utils::get_linear_index<memory_order>(nmp, local_mesh.dim);

  assert((nmp + Utils::Vector3i::broadcast(cao)) <= local_mesh.dim);
  for (int i = 0; i < cao; i++) {
    using Utils::bspline;

    ret.w_x[i] = bspline<cao>(i, dist[0]);
    ret.w_y[i] = bspline<cao>(i, dist[1]);
    ret.w_z[i] = bspline<cao>(i, dist[2]);
  }

  return ret;
}

/**
 * @brief P3M grid interpolation.
 *
 * This runs an kernel for every interpolation point
 * in a set of interpolation weights with the linear
 * grid index and the weight of the point as arguments.
 *
 * @param local_mesh Mesh info.
 * @param weights Set of weights
 * @param kernel The kernel to run.
 */
template <int cao, template <int> class WeightsStorage, class Kernel>
void p3m_interpolate(P3MLocalMesh const &local_mesh,
                     WeightsStorage<cao> const &weights, Kernel kernel) {
  auto q_ind = weights.ind;
  for (int i0 = 0; i0 < cao; i0++) {
    auto const w_x = weights.w_x[i0];
    for (int i1 = 0; i1 < cao; i1++) {
      auto const w_xy = w_x * weights.w_y[i1];
      for (int i2 = 0; i2 < cao; i2++) {
        auto const w_xyz = w_xy * weights.w_z[i2];
        kernel(q_ind, w_xyz);
        q_ind++;
      }
      q_ind += local_mesh.q_2_off;
    }
    q_ind += local_mesh.q_21_off;
  }
}
