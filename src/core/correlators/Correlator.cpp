/*
 Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

 This file is part of ESPResSo.

 ESPResSo is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ESPResSo is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Correlator.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

#include "utils/serialization/multi_array.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <limits>

namespace {
int min(int i, unsigned int j) { return std::min(i, static_cast<int>(j)); }
}

namespace Correlators {
/** The minimal version of compression function */
std::vector<double> compress_do_nothing(std::vector<double> const &A1,
                                        std::vector<double> const &A2) {
  return {};
}

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
                                   std::vector<double> const &B, Vector3d) {
  if (A.size() != B.size()) {
    throw std::runtime_error(
        "Error in scalar product: The vector sizes do not match");
  }

  return std::vector<double>(
      1, std::inner_product(A.begin(), A.end(), B.begin(), 0.0));
}

std::vector<double> componentwise_product(std::vector<double> const &A,
                                          std::vector<double> const &B,
                                          Vector3d) {
  std::vector<double> C(A.size());
  if (!(A.size() == B.size())) {
    throw std::runtime_error(
        "Error in componentwise product: The vector sizes do not match");
  }

  std::transform(A.begin(), A.end(), B.begin(), C.begin(),
                 std::multiplies<double>());

  return C;
}

std::vector<double> tensor_product(std::vector<double> const &A,
                                   std::vector<double> const &B, Vector3d) {
  std::vector<double> C(A.size() * B.size());
  auto C_it = C.begin();

  for (auto A_it = A.begin(); A_it != A.begin(); ++A_it) {
    for (auto B_it = B.begin(); B_it != B.end(); ++B_it) {
      *(C_it++) = *A_it * *B_it;
    }
  }
  assert(C_it == C.end());

  return C;
}

std::vector<double> square_distance_componentwise(std::vector<double> const &A,
                                                  std::vector<double> const &B,
                                                  Vector3d) {
  if (!(A.size() == B.size())) {
    throw std::runtime_error(
        "Error in square distance componentwise: The vector sizes do not "
        "match.");
  }

  std::vector<double> C(A.size());

  std::transform(
      A.begin(), A.end(), B.begin(), C.begin(),
      [](double a, double b) -> double { return (a - b) * (a - b); });

  return C;
}

// note: the argument name wsquare denotes that it value is w^2 while the user
// sets w
std::vector<double> fcs_acf(std::vector<double> const &A,
                            std::vector<double> const &B, Vector3d wsquare) {
  if (!(A.size() == B.size())) {
    throw std::runtime_error(
        "Error in square distance componentwise: The vector sizes do not "
        "match.");
  }

  auto const C_size = A.size() / 3;
  if (3 * C_size != A.size()) {
    throw std::runtime_error("Invalid dimensions.");
  }

  std::vector<double> C(C_size, 0);

  for (unsigned i = 0; i < C_size; i++) {
    for (int j = 0; j < 3; j++) {
      auto const &a = A[3 * i + j];
      auto const &b = B[3 * i + j];

      C[i] -= (a - b) * (a - b) / wsquare[j];
    }
  }

  std::transform(C.begin(), C.end(), C.begin(),
                 [](double c) -> double { return std::exp(c); });

  return C;
}

/* global variables */

/* Error codes */
constexpr const char init_errors[][64] = {
    "",                                                              // 0
    "No valid correlation given",                                    // 1
    "delta_t must be specified and > 0",                             // 2
    "tau_lin must be >= 2",                                          // 3
    "tau_max must be >= delta_t",                                    // 4
    "window_distance must be >1",                                    // 5
    "dimension of A was not >1",                                     // 6
    "dimension of B was not >1",                                     // 7
    "dimension of B must match dimension of A ",                     // 8
    "no proper function for first observable given",                 // 9
    "no proper function for second observable given",                // 10
    "no proper function for correlation operation given",            // 11
    "no proper function for compression of first observable given",  // 12
    "no proper function for compression of second observable given", // 13
    "tau_lin must be divisible by 2",                                // 14
    "dt is smaller than the MD timestep",                            // 15
    "dt is not a multiple of the MD timestep",                       // 16
    "cannot set compress2 for autocorrelation",                      // 17
    "diA must be divisible by 3 for fcs_acf",                        // 18
    "fcs_acf requires 3 additional parameters"                       // 19
};

int correlations_autoupdate = 0;

int Correlator::get_correlation_time(double *correlation_time) {
  // We calculate the correlation time for each m_dim_corr by normalizing the
  // correlation,
  // integrating it and finding out where C(tau)=tau;
  double C_tau;
  int ok_flag;
  for (unsigned j = 0; j < m_dim_corr; j++) {
    correlation_time[j] = 0.;
  }

  // here we still have to fix the stuff a bit!
  for (unsigned j = 0; j < m_dim_corr; j++) {
    C_tau = 1 * m_dt;
    ok_flag = 0;
    for (unsigned k = 1; k < m_n_result - 1; k++) {
      if (n_sweeps[k] == 0)
        break;
      C_tau += (result[k][j] / (double)n_sweeps[k] -
                A_accumulated_average[j] * B_accumulated_average[j] / n_data /
                    n_data) /
               (result[0][j] / n_sweeps[0]) * m_dt * (tau[k] - tau[k - 1]);

      if (exp(-tau[k] * m_dt / C_tau) + 2 * sqrt(tau[k] * m_dt / n_data) >
          exp(-tau[k - 1] * m_dt / C_tau) +
              2 * sqrt(tau[k - 1] * m_dt / n_data)) {
        correlation_time[j] =
            C_tau * (1 + (2 * (double)tau[k] + 1) / (double)n_data);
        ok_flag = 1;
        break;
      }
    }
    if (!ok_flag) {
      correlation_time[j] = -1;
    }
  }

  return 0;
}

void Correlator::initialize() {
  unsigned int i, j, k;
  hierarchy_depth = 0;
  // Class members are assigned via the initializer list

  // Input validation
  if (m_dt <= 0) {
    throw std::runtime_error(init_errors[2]);
  }

  if ((m_dt - time_step) < -1e-6 * time_step) {
    throw std::runtime_error(init_errors[15]);
  }

  // check if m_dt is a multiple of the md timestep
  if (std::abs(m_dt / time_step - round(m_dt / time_step)) > 1e-6) {
    throw std::runtime_error(init_errors[16]);
  }

  // Time steps and intervals
  update_frequency =
      std::floor(m_dt / time_step + std::numeric_limits<double>::round_error());

  if (m_tau_lin == 1) { // use the default
    m_tau_lin = (int)ceil(m_tau_max / m_dt);
    if (m_tau_lin % 2)
      m_tau_lin += 1;
  }

  if (m_tau_lin < 2) {
    throw std::runtime_error(init_errors[3]);
  }

  if (m_tau_lin % 2) {
    throw std::runtime_error(init_errors[14]);
  }

  if (m_tau_max <= m_dt) {
    throw std::runtime_error(init_errors[4]);

  } else { // set hierarchy depth which can  accommodate at least m_tau_max
    if ((m_tau_max / m_dt) < m_tau_lin) {
      hierarchy_depth = 1;
    } else {
      hierarchy_depth = (unsigned int)ceil(
          1 + log((m_tau_max / m_dt) / (m_tau_lin - 1)) / log(2.0));
    }
  }

  dim_A = 0;
  dim_B = 0;

  if (A_obs)
    dim_A = A_obs->n_values();
  if (!B_obs) {
    B_obs = A_obs;
  }

  dim_B = B_obs->n_values();

  if (dim_A < 1) {
    throw std::runtime_error(init_errors[6]);
  }

  // choose the correlation operation
  if (corr_operation_name == "") {
    throw std::runtime_error(init_errors[11]); // there is no reasonable default
  } else if (corr_operation_name == "componentwise_product") {
    m_dim_corr = dim_A;
    corr_operation = &componentwise_product;
    m_correlation_args = Vector3d{0, 0, 0};
  } else if (corr_operation_name == "tensor_product") {
    m_dim_corr = dim_A * dim_B;
    corr_operation = &tensor_product;
    m_correlation_args = Vector3d{0, 0, 0};
  } else if (corr_operation_name == "square_distance_componentwise") {
    m_dim_corr = dim_A;
    corr_operation = &square_distance_componentwise;
    m_correlation_args = Vector3d{0, 0, 0};
  } else if (corr_operation_name == "fcs_acf") {
    // note: user provides w=(wx,wy,wz) but we want to use
    // wsquare=(wx^2,wy^2,wz^2)
    if (m_correlation_args[0] <= 0 || m_correlation_args[1] <= 0 ||
        m_correlation_args[2] <= 0) {
      throw std::runtime_error("missing parameter for fcs_acf: w_x w_y w_z");
    }
    m_correlation_args[0] = m_correlation_args[0] * m_correlation_args[0];
    m_correlation_args[1] = m_correlation_args[1] * m_correlation_args[1];
    m_correlation_args[2] = m_correlation_args[2] * m_correlation_args[2];
    fprintf(stderr, "args2: %f %f %f\n", m_correlation_args[0],
            m_correlation_args[1], m_correlation_args[2]);
    if (dim_A % 3)
      throw std::runtime_error(init_errors[18]);
    m_dim_corr = dim_A / 3;
    corr_operation = &fcs_acf;
  } else if (corr_operation_name == "scalar_product") {
    m_dim_corr = 1;
    corr_operation = &scalar_product;
    m_correlation_args = Vector3d{0, 0, 0};
  } else {
    throw std::runtime_error(init_errors[11]);
  }

  // Choose the compression function
  if (compressA_name == "") { // this is the default
    compressA_name = "discard2";
    compressA = &compress_discard2;
  } else if (compressA_name == "discard2") {
    compressA = &compress_discard2;
  } else if (compressA_name == "discard1") {
    compressA = &compress_discard1;
  } else if (compressA_name == "linear") {
    compressA = &compress_linear;
  } else {
    throw std::runtime_error(init_errors[12]);
  }

  if (compressB_name == "") {
    compressB_name = compressA_name;
    compressB = compressA;
  } else if (compressB_name == "discard2") {
    compressB = &compress_discard2;
  } else if (compressB_name == "discard1") {
    compressB = &compress_discard1;
  } else if (compressB_name == "linear") {
    compressB = &compress_linear;
  } else {
    throw std::runtime_error(init_errors[13]);
  }

  A.resize(std::array<int, 2>{{hierarchy_depth, m_tau_lin + 1}});
  std::fill_n(A.data(), A.num_elements(), std::vector<double>(dim_A, 0));
  B.resize(std::array<int, 2>{{hierarchy_depth, m_tau_lin + 1}});
  std::fill_n(B.data(), B.num_elements(), std::vector<double>(dim_B, 0));

  n_data = 0;
  A_accumulated_average = std::vector<double>(dim_A, 0);
  B_accumulated_average = std::vector<double>(dim_B, 0);

  m_n_result = m_tau_lin + 1 + (m_tau_lin + 1) / 2 * (hierarchy_depth - 1);
  n_sweeps = std::vector<unsigned int>(m_n_result, 0);
  n_vals = std::vector<unsigned int>(hierarchy_depth, 0);

  result.resize(std::array<int, 2>{{m_n_result, m_dim_corr}});

  for (i = 0; i < m_n_result; i++) {
    for (j = 0; j < m_dim_corr; j++)
      // and initialize the values
      result[i][j] = 0;
  }

  newest = std::vector<unsigned int>(hierarchy_depth, m_tau_lin);

  tau.resize(m_n_result);
  for (i = 0; i < m_tau_lin + 1; i++) {
    tau[i] = i;
  }

  for (j = 1; j < hierarchy_depth; j++)
    for (k = 0; k < m_tau_lin / 2; k++) {
      tau[m_tau_lin + 1 + (j - 1) * m_tau_lin / 2 + k] =
          (k + (m_tau_lin / 2) + 1) * (1 << j);
    }
}

int Correlator::get_data() {
  if (finalized) {
    runtimeErrorMsg() << "No data can be added after finalize() was called.";
    return 0;
  }
  // We must now go through the hierarchy and make sure there is space for the
  // new
  // datapoint. For every hierarchy level we have to decide if it necessary to
  // move
  // something
  int i, j;
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;

  t++;

  highest_level_to_compress = -1;
  i = 0;
  j = 1;
  // Lets find out how far we have to go back in the hierarchy to make space for
  // the new value
  while (1) {
    if (((t - ((m_tau_lin + 1) * ((1 << (i + 1)) - 1) + 1)) % (1 << (i + 1)) ==
         0)) {
      if (i < (int(hierarchy_depth) - 1) && n_vals[i] > m_tau_lin) {
        highest_level_to_compress += 1;
        i++;
      } else
        break;
    } else
      break;
  }

  // Now we know we must make space on the levels 0..highest_level_to_compress
  // Now lets compress the data level by level.

  for (i = highest_level_to_compress; i >= 0; i--) {
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

  A[0][newest[0]] = A_obs->operator()(partCfg());
  if (A_obs != B_obs) {
    B[0][newest[0]] = B_obs->operator()(partCfg());
  } else {
    B[0][newest[0]] = A[0][newest[0]];
  }

  // Now we update the cumulated averages and variances of A and B
  n_data++;
  for (unsigned k = 0; k < dim_A; k++) {
    A_accumulated_average[k] += A[0][newest[0]][k];
  }

  for (unsigned k = 0; k < dim_B; k++) {
    B_accumulated_average[k] += B[0][newest[0]][k];
  }

  // Now update the lowest level correlation estimates
  for (j = 0; j < min(m_tau_lin + 1, n_vals[0]); j++) {
    index_new = newest[0];
    index_old = (newest[0] - j + m_tau_lin + 1) % (m_tau_lin + 1);
    auto const temp =
        (corr_operation)(A[0][index_old], B[0][index_new], m_correlation_args);
    assert(temp.size() == m_dim_corr);

    n_sweeps[j]++;
    for (unsigned k = 0; k < m_dim_corr; k++) {
      result[j][k] += temp[k];
    }
  }
  // Now for the higher ones
  for (int i = 1; i < highest_level_to_compress + 2; i++) {
    for (unsigned j = (m_tau_lin + 1) / 2 + 1;
         j < min(m_tau_lin + 1, n_vals[i]); j++) {
      index_new = newest[i];
      index_old = (newest[i] - j + m_tau_lin + 1) % (m_tau_lin + 1);
      index_res =
          m_tau_lin + (i - 1) * m_tau_lin / 2 + (j - m_tau_lin / 2 + 1) - 1;
      auto const temp = (corr_operation)(A[i][index_old], B[i][index_new],
                                         m_correlation_args);
      assert(temp.size() == m_dim_corr);

      n_sweeps[index_res]++;
      for (unsigned k = 0; k < m_dim_corr; k++) {
        result[index_res][k] += temp[k];
      }
    }
  }

  m_last_update = sim_time;
  return 0;
}

int Correlator::finalize() {
  if (finalized) {
    runtimeErrorMsg() << "Correlator::finalize() can only be called once.";
    return 0;
  }
  // We must now go through the hierarchy and make sure there is space for the
  // new
  // datapoint. For every hierarchy level we have to decide if it necessary to
  // move
  // something
  int i, j;
  int ll = 0;      // current lowest level
  int vals_ll = 0; // number of values remaining in the lowest level
  int highest_level_to_compress;
  unsigned int index_new, index_old, index_res;

  // make a flag that the correlation is finalized
  finalized = 1;

  for (ll = 0; ll < hierarchy_depth - 1; ll++) {
    if (n_vals[ll] > m_tau_lin + 1)
      vals_ll = m_tau_lin + n_vals[ll] % 2;
    else
      vals_ll = n_vals[ll];

    while (vals_ll) {
      // Check, if we will want to push the value from the lowest level
      if (vals_ll % 2) {
        highest_level_to_compress = ll;
      } else {
        highest_level_to_compress = -1;
      }

      i = ll + 1; // lowest level, for which we have to check for compression
      j = 1;
      // Lets find out how far we have to go back in the hierarchy to make space
      // for the new value
      while (highest_level_to_compress > -1) {
        // printf("test level %d for compression, n_vals=%d ... ",i,n_vals[i]);
        if (n_vals[i] % 2) {
          if (i < (hierarchy_depth - 1) && n_vals[i] > m_tau_lin) {
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
      // Now lets compress the data level by level.

      for (i = highest_level_to_compress; i >= ll; i--) {
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
      for (i = ll + 1; i < highest_level_to_compress + 2; i++) {
        for (j = (m_tau_lin + 1) / 2 + 1;
             j < min(m_tau_lin + 1, n_vals[i]); j++) {
          index_new = newest[i];
          index_old = (newest[i] - j + m_tau_lin + 1) % (m_tau_lin + 1);
          index_res =
              m_tau_lin + (i - 1) * m_tau_lin / 2 + (j - m_tau_lin / 2 + 1) - 1;

          auto const temp = (corr_operation)(A[i][index_old], B[i][index_new],
                                             m_correlation_args);
          assert(temp.size() == m_dim_corr);

          n_sweeps[index_res]++;
          for (unsigned k = 0; k < m_dim_corr; k++) {
            result[index_res][k] += temp[k];
          }
        }
      }
    }
  }
  return 0;
}

void Correlator::start_auto_update() {
  if (update_frequency > 0) {
    correlations_autoupdate = 1;
    autoupdate = 1;
    m_last_update = sim_time;
  } else {
    throw std::runtime_error(
        "Could not start autoupdate: update frequency not set");
  }
}

void Correlator::stop_auto_update() {
  autoupdate = 0;
  // Todo
  // Insert logic to determine if global correlations_auto_update can be set to
  // 0
}

std::vector<double> Correlator::get_correlation() {
  std::vector<double> res;

  // time + n_sweeps + corr_1...corr_n
  int cols = 2 + m_dim_corr;
  res.resize(m_n_result * cols);

  for (int i = 0; i < m_n_result; i++) {
    res[cols * i + 0] = tau[i] * m_dt;
    res[cols * i + 1] = n_sweeps[i];
    for (int k = 0; k < m_dim_corr; k++) {
      if (n_sweeps[i] > 0) {
        res[cols * i + 2 + k] = result[i][k] / n_sweeps[i];
      } else {
        res[cols * i + 2 + k] = 0;
      }
    }
  }
  return res;
}

std::string Correlator::get_internal_state() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << t;
  oa << m_n_result;
  oa << A;
  oa << B;
  oa << result;
  oa << n_sweeps;
  oa << n_vals;
  oa << newest;
  oa << A_accumulated_average;
  oa << B_accumulated_average;
  oa << n_data;
  oa << m_last_update;

  return ss.str();
}
void Correlator::set_internal_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  ia >> t;
  ia >> m_n_result;
  ia >> A;
  ia >> B;
  ia >> result;
  ia >> n_sweeps;
  ia >> n_vals;
  ia >> newest;
  ia >> A_accumulated_average;
  ia >> B_accumulated_average;
  ia >> n_data;
  ia >> m_last_update;
}

} // Namespace Correlators
