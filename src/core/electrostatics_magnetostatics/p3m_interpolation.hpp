/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef ESPRESSO_P3M_INTERPOLATION_HPP
#define ESPRESSO_P3M_INTERPOLATION_HPP

#include <utils/Span.hpp>
#include <utils/index.hpp>
#include <utils/math/bspline.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <cassert>
#include <tuple>
#include <vector>

/**
 * @brief Interpolation weights for one point.
 *
 * Interpolation weights and  grid offset for one point.
 *
 * @tparam cao Interpolation order.
 */
template <int cao> struct InterpolationWeights {
  /** Linear index of the corner of the interpolation cube. */
  int ind;
  /** Weights for the directions */
  Utils::Array<double, cao> w_x, w_y, w_z;
};

/**
 * @brief Cache for interpolation weights.
 *
 * This is a storage container for interpolation weights of
 * type InterpolationWeights.
 */
class p3m_interpolation_cache {
  size_t m_cao = 0;
  /** Charge fractions for mesh assignment. */
  std::vector<double> ca_frac;
  /** index of first mesh point for charge assignment. */
  std::vector<int> ca_fmp;

public:
  /**
   * @brief Number of points in the cache.
   * @return Number of points currently in the cache.
   */
  auto size() const { return ca_fmp.size(); }

  /**
   * @brief Charge assignment order the weights are for.
   * @return The charge assignment order.
   */
  auto cao() const { return m_cao; }

  /**
   * @brief Push back weights for one point.
   *
   * @tparam cao Interpolation order has to match the order
   *         set at last call to @ref p3m_interpolation_cache::reset.
   * @param w Interpolation weights to store.
   */
  template <int cao> void store(const InterpolationWeights<cao> &w) {
    assert(cao == m_cao);

    ca_fmp.push_back(w.ind);
    auto it = std::back_inserter(ca_frac);
    boost::copy(w.w_x, it);
    boost::copy(w.w_y, it);
    boost::copy(w.w_z, it);
  }

  /**
   * @brief Load entry from the cache.
   *
   * This loads an entry at an index from the cache,
   * the entries are indexed by the order they were stored.
   *
   * @tparam cao Interpolation order has to match the order
   *         set at last call to @ref p3m_interpolation_cache::reset.
   * @param i Index of the entry to load.
   * @return i-it interpolation weights.
   */
  template <int cao> InterpolationWeights<cao> load(size_t i) const {
    assert(cao == m_cao);

    using Utils::make_const_span;
    assert(i < size());

    InterpolationWeights<cao> ret;
    ret.ind = ca_fmp[i];

    auto const offset = ca_frac.data() + 3 * i * m_cao;
    boost::copy(make_const_span(offset + 0 * m_cao, m_cao), ret.w_x.begin());
    boost::copy(make_const_span(offset + 1 * m_cao, m_cao), ret.w_y.begin());
    boost::copy(make_const_span(offset + 2 * m_cao, m_cao), ret.w_z.begin());

    return ret;
  }

  /**
   * @brief Reset the cache.
   *
   * @param cao Interpolation order.
   */
  void reset(int cao) {
    m_cao = cao;
    ca_frac.clear();
    ca_fmp.clear();
  }
};

/**
 * @brief Calculate the P-th order interpolation weights.
 *
 * As described in from @cite hockney88a 5-189 (or 8-61).
 * The weights are also tabulated in @cite deserno98a @cite deserno98b.
 */
template <int cao>
InterpolationWeights<cao>
p3m_calculate_interpolation_weights(const Utils::Vector3d &position,
                                    const Utils::Vector3d &ai,
                                    p3m_local_mesh const &local_mesh) {
  /** position shift for calc. of first assignment mesh point. */
  static auto const pos_shift = std::floor((cao - 1) / 2.0) - (cao % 2) / 2.0;

  /* distance to nearest mesh point */
  Utils::Vector3d dist;

  /* nearest mesh point */
  Utils::Vector3i nmp;

  for (int d = 0; d < 3; d++) {
    /* particle position in mesh coordinates */
    auto const pos = ((position[d] - local_mesh.ld_pos[d]) * ai[d]) - pos_shift;

    nmp[d] = (int)pos;

    /* distance to nearest mesh point */
    dist[d] = (pos - nmp[d]) - 0.5;
  }

  InterpolationWeights<cao> ret;

  /* 3d-array index of nearest mesh point */
  ret.ind = Utils::get_linear_index(nmp, local_mesh.dim,
                                    Utils::MemoryOrder::ROW_MAJOR);
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
template <int cao, class Kernel>
void p3m_interpolate(p3m_local_mesh const &local_mesh,
                     InterpolationWeights<cao> const &weights, Kernel kernel) {
  auto q_ind = weights.ind;
  for (int i0 = 0; i0 < cao; i0++) {
    auto const tmp0 = weights.w_x[i0];
    for (int i1 = 0; i1 < cao; i1++) {
      auto const tmp1 = tmp0 * weights.w_y[i1];
      for (int i2 = 0; i2 < cao; i2++) {
        kernel(q_ind, tmp1 * weights.w_z[i2]);

        q_ind++;
      }
      q_ind += local_mesh.q_2_off;
    }
    q_ind += local_mesh.q_21_off;
  }
}

#endif // ESPRESSO_P3M_INTERPOLATION_HPP
