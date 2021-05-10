/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "Correlator.hpp"

#include "integrate.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>
#include <utils/serialization/multi_array.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace {
int min(int i, unsigned int j) { return std::min(i, static_cast<int>(j)); }
} // namespace

namespace Accumulators {
/** Compress computing arithmetic mean: A_compressed=(A1+A2)/2 */
std::vector<double> compress_linear(std::vector<double> const &A1,
                                    std::vector<double> const &A2) {
  assert(A1.size() == A2.size());
  std::vector<double> A_compressed(A1.size());

  std::transform(A1.begin(), A1.end(), A2.begin(), A_compressed.begin(),
                 [](double a, double b) -> double { return 0.5 * (a + b); });

  return A_compressed;
}

/** Compress discarding the 1st argument and return the 2nd */
std::vector<double> compress_discard1(std::vector<double> const &A1,
                                      std::vector<double> const &A2) {
  assert(A1.size() == A2.size());
  std::vector<double> A_compressed(A2);
  return A_compressed;
}

/** Compress discarding the 2nd argument and return the 1st */
std::vector<double> compress_discard2(std::vector<double> const &A1,
                                      std::vector<double> const &A2) {
  assert(A1.size() == A2.size());
  std::vector<double> A_compressed(A1);
  return A_compressed;
}

std::vector<double> scalar_product(std::vector<double> const &A,
                                   std::vector<double> const &B,
                                   Utils::Vector3d const &) {
  if (A.size() != B.size()) {
    throw std::runtime_error(
        "Error in scalar product: The vector sizes do not match");
  }

  return std::vector<double>(
      1, std::inner_product(A.begin(), A.end(), B.begin(), 0.0));
}

std::vector<double> componentwise_product(std::vector<double> const &A,
                                          std::vector<double> const &B,
                                          Utils::Vector3d const &) {
  std::vector<double> C(A.size());
  if (A.size() != B.size()) {
    throw std::runtime_error(
        "Error in componentwise product: The vector sizes do not match");
  }

  std::transform(A.begin(), A.end(), B.begin(), C.begin(), std::multiplies<>());

  return C;
}

std::vector<double> tensor_product(std::vector<double> const &A,
                                   std::vector<double> const &B,
                                   Utils::Vector3d const &) {
  std::vector<double> C(A.size() * B.size());
  auto C_it = C.begin();

  for (double a : A) {
    for (double b : B) {
      *(C_it++) = a * b;
    }
  }

  return C;
}

std::vector<double> square_distance_componentwise(std::vector<double> const &A,
                                                  std::vector<double> const &B,
                                                  Utils::Vector3d const &) {
  if (A.size() != B.size()) {
    throw std::runtime_error(
        "Error in square distance componentwise: The vector sizes do not "
        "match.");
  }

  std::vector<double> C(A.size());

  std::transform(
      A.begin(), A.end(), B.begin(), C.begin(),
      [](double a, double b) -> double { return Utils::sqr(a - b); });

  return C;
}

// note: the argument name wsquare denotes that its value is w^2 while the user
// sets w
std::vector<double> fcs_acf(std::vector<double> const &A,
                            std::vector<double> const &B,
                            Utils::Vector3d const &wsquare) {
  if (A.size() != B.size()) {
    throw std::runtime_error(
        "Error in fcs_acf: The vector sizes do not match.");
  }

  auto const C_size = A.size() / 3;
  assert(3 * C_size == A.size());

  std::vector<double> C(C_size, 0);

  for (size_t i = 0; i < C_size; i++) {
    for (int j = 0; j < 3; j++) {
      auto const &a = A[3 * i + j];
      auto const &b = B[3 * i + j];

      C[i] -= Utils::sqr(a - b) / wsquare[j];
    }
  }

  std::transform(C.begin(), C.end(), C.begin(),
                 [](double c) -> double { return std::exp(c); });

  return C;
}

void Correlator::initialize() {
  // Class members are assigned via the initializer list

  if (m_tau_lin == 1) { // use the default
    m_tau_lin = static_cast<int>(ceil(m_tau_max / m_dt));
    if (m_tau_lin % 2)
      m_tau_lin += 1;
  }

  if (m_tau_lin < 2) {
    throw std::runtime_error("tau_lin must be >= 2");
  }

  if (m_tau_lin % 2) {
    throw std::runtime_error("tau_lin must be divisible by 2");
  }

  if (m_tau_max <= m_dt) {
    throw std::runtime_error("tau_max must be >= delta_t (delta_N too large)");
  }
  // set hierarchy depth which can accommodate at least m_tau_max
  if ((m_tau_max / m_dt) < m_tau_lin) {
    m_hierarchy_depth = 1;
  } else {
    m_hierarchy_depth = static_cast<int>(
        ceil(1 + log((m_tau_max / m_dt) / (m_tau_lin - 1)) / log(2.0)));
  }

  assert(A_obs);
  assert(B_obs);
  dim_A = A_obs->n_values();
  dim_B = B_obs->n_values();

  if (dim_A == 0) {
    throw std::runtime_error("dimension of first observable has to be >= 1");
  }
  if (dim_B == 0) {
    throw std::runtime_error("dimension of second observable has to be >= 1");
  }

  // choose the correlation operation
  if (corr_operation_name == "componentwise_product") {
    m_dim_corr = dim_A;
    m_shape = A_obs->shape();
    corr_operation = &componentwise_product;
    m_correlation_args = Utils::Vector3d{0, 0, 0};
  } else if (corr_operation_name == "tensor_product") {
    m_dim_corr = dim_A * dim_B;
    m_shape = {dim_A, dim_B};
    corr_operation = &tensor_product;
    m_correlation_args = Utils::Vector3d{0, 0, 0};
  } else if (corr_operation_name == "square_distance_componentwise") {
    m_dim_corr = dim_A;
    m_shape = A_obs->shape();
    corr_operation = &square_distance_componentwise;
    m_correlation_args = Utils::Vector3d{0, 0, 0};
  } else if (corr_operation_name == "fcs_acf") {
    // note: user provides w=(wx,wy,wz) but we want to use
    // wsquare=(wx^2,wy^2,wz^2)
    if (m_correlation_args[0] <= 0 || m_correlation_args[1] <= 0 ||
        m_correlation_args[2] <= 0) {
      throw std::runtime_error("missing parameter for fcs_acf: w_x w_y w_z");
    }
    m_correlation_args =
        Utils::hadamard_product(m_correlation_args, m_correlation_args);
    if (dim_A % 3)
      throw std::runtime_error("dimA must be divisible by 3 for fcs_acf");
    m_dim_corr = dim_A / 3;
    m_shape = A_obs->shape();
    if (m_shape.back() != 3)
      throw std::runtime_error(
          "the last dimension of dimA must be 3 for fcs_acf");
    m_shape.pop_back();
    corr_operation = &fcs_acf;
  } else if (corr_operation_name == "scalar_product") {
    m_dim_corr = 1;
    m_shape = {1};
    corr_operation = &scalar_product;
    m_correlation_args = Utils::Vector3d{0, 0, 0};
  } else {
    throw std::invalid_argument("correlation operation '" +
                                corr_operation_name + "' not implemented");
  }

  // Choose the compression function
  if (compressA_name == "discard2") {
    compressA = &compress_discard2;
  } else if (compressA_name == "discard1") {
    compressA = &compress_discard1;
  } else if (compressA_name == "linear") {
    compressA = &compress_linear;
  } else {
    throw std::invalid_argument("unknown compression method '" +
                                compressA_name + "' for first observable");
  }

  if (compressB_name == "discard2") {
    compressB = &compress_discard2;
  } else if (compressB_name == "discard1") {
    compressB = &compress_discard1;
  } else if (compressB_name == "linear") {
    compressB = &compress_linear;
  } else {
    throw std::invalid_argument("unknown compression method '" +
                                compressB_name + "' for second observable");
  }

  A.resize(std::array<int, 2>{{m_hierarchy_depth, m_tau_lin + 1}});
  std::fill_n(A.data(), A.num_elements(), std::vector<double>(dim_A, 0));
  B.resize(std::array<int, 2>{{m_hierarchy_depth, m_tau_lin + 1}});
  std::fill_n(B.data(), B.num_elements(), std::vector<double>(dim_B, 0));

  n_data = 0;
  A_accumulated_average = std::vector<double>(dim_A, 0);
  B_accumulated_average = std::vector<double>(dim_B, 0);

  auto const n_result = n_values();
  n_sweeps = std::vector<size_t>(n_result, 0);
  n_vals = std::vector<unsigned int>(m_hierarchy_depth, 0);

  result.resize(std::array<size_t, 2>{{n_result, m_dim_corr}});

  for (size_t i = 0; i < n_result; i++) {
    for (size_t j = 0; j < m_dim_corr; j++) {
      // and initialize the values
      result[i][j] = 0;
    }
  }

  newest = std::vector<size_t>(m_hierarchy_depth, m_tau_lin);

  tau.resize(n_result);
  for (int i = 0; i < m_tau_lin + 1; i++) {
    tau[i] = i;
  }

  for (int j = 1; j < m_hierarchy_depth; j++) {
    for (int k = 0; k < m_tau_lin / 2; k++) {
      tau[m_tau_lin + 1 + (j - 1) * m_tau_lin / 2 + k] =
          (k + (m_tau_lin / 2) + 1) * (1 << j);
    }
  }
}

void Correlator::update() {
  if (finalized) {
    throw std::runtime_error(
        "No data can be added after finalize() was called.");
  }
  // We must now go through the hierarchy and make sure there is space for the
  // new datapoint. For every hierarchy level we have to decide if it is
  // necessary to move something
  int highest_level_to_compress = -1;

  t++;

  // Let's find out how far we have to go back in the hierarchy to make space
  // for the new value
  int i = 0;
  while (true) {
    if (((t - ((m_tau_lin + 1) * ((1 << (i + 1)) - 1) + 1)) % (1 << (i + 1)) ==
         0)) {
      if (i < (m_hierarchy_depth - 1) && n_vals[i] > m_tau_lin) {
        highest_level_to_compress += 1;
        i++;
      } else
        break;
    } else
      break;
  }

  // Now we know we must make space on the levels 0..highest_level_to_compress
  // Now let's compress the data level by level.

  for (int i = highest_level_to_compress; i >= 0; i--) {
    // We increase the index indicating the newest on level i+1 by one (plus
    // folding)
    newest[i + 1] = (newest[i + 1] + 1) % (m_tau_lin + 1);
    n_vals[i + 1] += 1;
    A[i + 1][newest[i + 1]] =
        (*compressA)(A[i][(newest[i] + 1) % (m_tau_lin + 1)],
                     A[i][(newest[i] + 2) % (m_tau_lin + 1)]);
    B[i + 1][newest[i + 1]] =
        (*compressB)(B[i][(newest[i] + 1) % (m_tau_lin + 1)],
                     B[i][(newest[i] + 2) % (m_tau_lin + 1)]);
  }

  newest[0] = (newest[0] + 1) % (m_tau_lin + 1);
  n_vals[0]++;

  A[0][newest[0]] = A_obs->operator()();
  if (A_obs != B_obs) {
    B[0][newest[0]] = B_obs->operator()();
  } else {
    B[0][newest[0]] = A[0][newest[0]];
  }

  // Now we update the cumulated averages and variances of A and B
  n_data++;
  for (size_t k = 0; k < dim_A; k++) {
    A_accumulated_average[k] += A[0][newest[0]][k];
  }

  for (size_t k = 0; k < dim_B; k++) {
    B_accumulated_average[k] += B[0][newest[0]][k];
  }

  // Now update the lowest level correlation estimates
  for (unsigned j = 0; j < min(m_tau_lin + 1, n_vals[0]); j++) {
    auto const index_new = newest[0];
    auto const index_old = (newest[0] - j + m_tau_lin + 1) % (m_tau_lin + 1);
    auto const temp =
        (corr_operation)(A[0][index_old], B[0][index_new], m_correlation_args);
    assert(temp.size() == m_dim_corr);

    n_sweeps[j]++;
    for (size_t k = 0; k < m_dim_corr; k++) {
      result[j][k] += temp[k];
    }
  }
  // Now for the higher ones
  for (int i = 1; i < highest_level_to_compress + 2; i++) {
    for (unsigned j = (m_tau_lin + 1) / 2 + 1;
         j < min(m_tau_lin + 1, n_vals[i]); j++) {
      auto const index_new = newest[i];
      auto const index_old = (newest[i] - j + m_tau_lin + 1) % (m_tau_lin + 1);
      auto const index_res =
          m_tau_lin + (i - 1) * m_tau_lin / 2 + (j - m_tau_lin / 2 + 1) - 1;
      auto const temp = (corr_operation)(A[i][index_old], B[i][index_new],
                                         m_correlation_args);
      assert(temp.size() == m_dim_corr);

      n_sweeps[index_res]++;
      for (size_t k = 0; k < m_dim_corr; k++) {
        result[index_res][k] += temp[k];
      }
    }
  }
}

int Correlator::finalize() {
  if (finalized) {
    throw std::runtime_error("Correlator::finalize() can only be called once.");
  }
  // We must now go through the hierarchy and make sure there is space for the
  // new datapoint. For every hierarchy level we have to decide if it is
  // necessary to move something

  // mark the correlation as finalized
  finalized = true;

  for (int ll = 0; ll < m_hierarchy_depth - 1; ll++) {
    int vals_ll; // number of values remaining in the lowest level
    if (n_vals[ll] > m_tau_lin + 1)
      vals_ll = m_tau_lin + static_cast<int>(n_vals[ll]) % 2;
    else
      vals_ll = n_vals[ll];

    while (vals_ll) {
      // Check, if we will want to push the value from the lowest level
      int highest_level_to_compress = -1;
      if (vals_ll % 2) {
        highest_level_to_compress = ll;
      }

      int i = ll + 1; // lowest level for which we have to check for compression
      // Let's find out how far we have to go back in the hierarchy to make
      // space for the new value
      while (highest_level_to_compress > -1) {
        if (n_vals[i] % 2) {
          if (i < (m_hierarchy_depth - 1) && n_vals[i] > m_tau_lin) {
            highest_level_to_compress += 1;
            i++;
          } else {
            break;
          }
        } else {
          break;
        }
      }
      vals_ll -= 1;

      // Now we know we must make space on the levels
      // 0..highest_level_to_compress
      // Now let's compress the data level by level.

      for (int i = highest_level_to_compress; i >= ll; i--) {
        // We increase the index indicating the newest on level i+1 by one (plus
        // folding)
        newest[i + 1] = (newest[i + 1] + 1) % (m_tau_lin + 1);
        n_vals[i + 1] += 1;

        (*compressA)(A[i][(newest[i] + 1) % (m_tau_lin + 1)],
                     A[i][(newest[i] + 2) % (m_tau_lin + 1)]);
        (*compressB)(B[i][(newest[i] + 1) % (m_tau_lin + 1)],
                     B[i][(newest[i] + 2) % (m_tau_lin + 1)]);
      }
      newest[ll] = (newest[ll] + 1) % (m_tau_lin + 1);

      // We only need to update correlation estimates for the higher levels
      for (int i = ll + 1; i < highest_level_to_compress + 2; i++) {
        for (int j = (m_tau_lin + 1) / 2 + 1; j < min(m_tau_lin + 1, n_vals[i]);
             j++) {
          auto const index_new = newest[i];
          auto const index_old =
              (newest[i] - j + m_tau_lin + 1) % (m_tau_lin + 1);
          auto const index_res =
              m_tau_lin + (i - 1) * m_tau_lin / 2 + (j - m_tau_lin / 2 + 1) - 1;

          auto const temp = (corr_operation)(A[i][index_old], B[i][index_new],
                                             m_correlation_args);
          assert(temp.size() == m_dim_corr);

          n_sweeps[index_res]++;
          for (size_t k = 0; k < m_dim_corr; k++) {
            result[index_res][k] += temp[k];
          }
        }
      }
    }
  }
  return 0;
}

std::vector<double> Correlator::get_correlation() {
  auto const n_result = n_values();
  std::vector<double> res(n_result * m_dim_corr);

  for (size_t i = 0; i < n_result; i++) {
    auto const index = m_dim_corr * i;
    for (size_t k = 0; k < m_dim_corr; k++) {
      if (n_sweeps[i]) {
        res[index + k] = result[i][k] / static_cast<double>(n_sweeps[i]);
      }
    }
  }
  return res;
}

std::vector<double> Correlator::get_lag_times() const {
  std::vector<double> res(n_values());
  boost::transform(tau, res.begin(),
                   [dt = m_dt](auto const &a) { return a * dt; });
  return res;
}

std::string Correlator::get_internal_state() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << t;
  oa << m_shape;
  oa << A;
  oa << B;
  oa << result;
  oa << n_sweeps;
  oa << n_vals;
  oa << newest;
  oa << A_accumulated_average;
  oa << B_accumulated_average;
  oa << n_data;

  return ss.str();
}

void Correlator::set_internal_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  ia >> t;
  ia >> m_shape;
  ia >> A;
  ia >> B;
  ia >> result;
  ia >> n_sweeps;
  ia >> n_vals;
  ia >> newest;
  ia >> A_accumulated_average;
  ia >> B_accumulated_average;
  ia >> n_data;
}

} // namespace Accumulators
