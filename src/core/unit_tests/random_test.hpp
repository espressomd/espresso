/*
 * Copyright (C) 2019 The ESPResSo project
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
#ifndef CORE_UNIT_TESTS_RANDOM_TEST_HPP
#define CORE_UNIT_TESTS_RANDOM_TEST_HPP

/* Helper functions to compute random numbers covariance in a single pass */

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variates/covariate.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/variant.hpp>
#include <functional>
#include <tuple>

#include <utils/Vector.hpp>

#include "random.hpp"

namespace Utils {
using VariantVectorXd = boost::variant<double, Vector2d, Vector3d, Vector4d>;
} // namespace Utils

using Utils::VariantVectorXd;

namespace {

using Utils::Vector;

class visitor_size : public boost::static_visitor<size_t> {
public:
  template <size_t N> size_t operator()(Vector<double, N> const &v) const {
    return v.size();
  }
  size_t operator()(double v) const { return 1; }
};

class visitor_get : public boost::static_visitor<double> {
public:
  template <size_t N>
  double operator()(Vector<double, N> const &v, size_t i) const {
    return v[i];
  }
  double operator()(double v, size_t i) const {
    assert(i == 0);
    return v;
  }
};

size_t get_size(VariantVectorXd const &vec) {
  return boost::apply_visitor(visitor_size(), vec);
}

double get_value(VariantVectorXd const &vec, size_t i) {
  return boost::apply_visitor(
      std::bind(visitor_get(), std::placeholders::_1, i), vec);
}

template <typename T> auto square_matrix(size_t N) {
  return std::vector<std::vector<T>>(N, std::vector<T>(N));
}

} // namespace

/** Draw a large sample of 3D vectors from PRNGs and compute statistics.
 *  Parameter @p noise_function is a generator that returns @f$ N @f$ vectors
 *  of size @f$ M_i @f$. The following statistics are evaluated: @f$ N @f$ means
 *  and @f$ N @f$ variances (samples are uncorrelated across axes, so pooling
 *  them is fine), and a covariance and a correlation matrix of size
 *  @f$ \sum M_i @f$.
 */
std::tuple<std::vector<double>, std::vector<double>,
           std::vector<std::vector<double>>, std::vector<std::vector<double>>>
noise_statistics(std::function<std::vector<VariantVectorXd>()> noise_function,
                 size_t sample_size) {

  // get size of the arrays and size of the triangular correlation matrix
  auto const first_value = noise_function();
  auto const n_vectors = first_value.size();
  std::vector<size_t> dimensions(n_vectors);
  std::transform(first_value.begin(), first_value.end(), dimensions.begin(),
                 [](auto const &element) { return get_size(element); });
  auto const matrix_dim = std::accumulate(dimensions.begin(), dimensions.end(),
                                          0, std::plus<size_t>());

  // set up boost accumulators
  namespace ba = boost::accumulators;
  namespace bt = boost::accumulators::tag;
  using stat_variance = ba::stats<bt::mean, bt::variance(ba::lazy)>;
  using stat_covariance = ba::stats<bt::covariance<double, bt::covariate1>>;
  using boost_variance = ba::accumulator_set<double, stat_variance>;
  using boost_covariance = ba::accumulator_set<double, stat_covariance>;
  std::vector<boost_variance> acc_variance(n_vectors);
  auto acc_covariance = ::square_matrix<boost_covariance>(matrix_dim);

  // accumulate
  for (size_t step = 0; step < sample_size; ++step) {
    auto const noise_tuple = noise_function();
    // for each vector, pool the random numbers of all columns
    for (size_t vec1 = 0; vec1 < dimensions.size(); ++vec1) {
      for (size_t col1 = 0; col1 < dimensions[vec1]; ++col1) {
        acc_variance[vec1](::get_value(noise_tuple[vec1], col1));
      }
    }
    // fill the covariance matrix (upper triangle)
    size_t index1 = 0;
    for (size_t vec1 = 0; vec1 < dimensions.size(); ++vec1) {
      for (size_t col1 = 0; col1 < dimensions[vec1]; ++col1) {
        size_t index2 = index1;
        for (size_t vec2 = vec1; vec2 < dimensions.size(); ++vec2) {
          for (size_t col2 = (vec2 == vec1) ? col1 : 0; col2 < dimensions[vec2];
               ++col2) {
            acc_covariance[index1][index2](
                ::get_value(noise_tuple[vec1], col1),
                ba::covariate1 = ::get_value(noise_tuple[vec2], col2));
            index2++;
          }
        }
        index1++;
      }
    }
  }

  // compute statistics
  std::vector<double> means(n_vectors);
  std::vector<double> variances(n_vectors);
  for (size_t i = 0; i < n_vectors; ++i) {
    means[i] = ba::mean(acc_variance[i]);
    variances[i] = ba::variance(acc_variance[i]);
  }
  auto covariance = ::square_matrix<double>(matrix_dim);
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = i; j < matrix_dim; ++j) {
      covariance[i][j] = covariance[j][i] =
          ba::covariance(acc_covariance[i][j]);
    }
  }
  auto correlation = ::square_matrix<double>(matrix_dim);
  for (size_t i = 0; i < matrix_dim; ++i) {
    for (size_t j = i; j < matrix_dim; ++j) {
      correlation[i][j] = correlation[j][i] =
          covariance[i][j] / sqrt(covariance[i][i] * covariance[j][j]);
    }
  }

  return std::make_tuple(means, variances, covariance, correlation);
}

boost::test_tools::predicate_result correlation_almost_equal(
    std::vector<std::vector<double>> const &correlation_matrix, size_t i,
    size_t j, double reference, double threshold) {
  auto const value = correlation_matrix[i][j];
  auto const diff = std::abs(value - reference);
  if (diff > threshold) {
    boost::test_tools::predicate_result res(false);
    res.message() << "The correlation coefficient M[" << i << "][" << j << "]{"
                  << value << "} differs from " << reference << " by " << diff
                  << " (> " << threshold << ")";
    return res;
  }
  return true;
}
#endif
