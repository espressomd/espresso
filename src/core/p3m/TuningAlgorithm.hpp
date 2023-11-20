/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "config/config.hpp"

#if defined(P3M) || defined(DP3M)

#include "p3m/TuningLogger.hpp"
#include "p3m/common.hpp"

#include <utils/Vector.hpp>

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <tuple>

namespace System {
class System;
}

/**
 * @brief Tuning algorithm for P3M.
 *
 * The algorithm basically determines the mesh, cao
 * and then the real-space cutoff, in this order.
 *
 * For each mesh, the optimal cao for the previous mesh is re-used as an
 * initial guess, and the algorithm checks whether increasing or decreasing
 * it leads to a better solution. This is efficient, since the optimal cao
 * only changes little with the meshes in general.
 *
 * The real-space cutoff for a given mesh and cao is determined via a
 * bisection on the error estimate, which determines where the error
 * estimate equals the required accuracy. Therefore the smallest possible,
 * i.e. fastest real-space cutoff is determined.
 *
 * Both the search over mesh and cao stop to search in a specific
 * direction once the computation time is significantly higher
 * than the currently known optimum.
 */
class TuningAlgorithm {
protected:
  System::System &m_system;

private:
  int m_timings;
  std::size_t m_n_trials;

protected:
  double m_prefactor;
  std::unique_ptr<TuningLogger> m_logger = nullptr;
  int cao_min = -1, cao_max = -1, cao_best = -1;
  double m_r_cut_iL_min = -1., m_r_cut_iL_max = -1.;

  /**
   * @brief Granularity of the time measurement (milliseconds).
   * Tuning halts when the runtime is larger than the best time plus this value.
   */
  static auto constexpr time_granularity = 2.;

  /**
   * @brief Maximal number of consecutive trials that don't improve runtime.
   * Tuning halts when this threshold is reached.
   */
  static auto constexpr max_n_consecutive_trials = 20;

  /** @brief Value for invalid time measurements. */
  static auto constexpr time_sentinel = std::numeric_limits<double>::max();

public:
  TuningAlgorithm(System::System &system, double prefactor, int timings)
      : m_system{system}, m_timings{timings}, m_n_trials{0ul}, m_prefactor{
                                                                   prefactor} {}

  virtual ~TuningAlgorithm() = default;

  struct Parameters {
    Utils::Vector3i mesh = {};
    int cao = -1;
    double alpha_L = -1.;
    double r_cut_iL = -1.;
    double accuracy = -1.;
    double time = std::numeric_limits<double>::max();
  };

  /** @brief Get the P3M parameters. */
  virtual P3MParameters &get_params() = 0;

  /** @brief Re-initialize the currently active solver. */
  virtual void on_solver_change() const = 0;

  /** @brief Tuning loop entry point. */
  virtual TuningAlgorithm::Parameters get_time() = 0;

  /** @brief Configure the logger. */
  virtual void setup_logger(bool verbose) = 0;

  /** @brief Determine a sensible range for the mesh. */
  virtual void determine_mesh_limits() = 0;

  /** @brief Determine a sensible range for the real-space cutoff. */
  void determine_r_cut_limits();

  /** @brief Determine a sensible range for the charge assignment order. */
  void determine_cao_limits(int initial_cao);

  /**
   * @brief Get the minimal error for this combination of parameters.
   *
   * The real-space error is tuned such that it contributes half of the
   * total error, and then the k-space error is calculated.
   * If an optimal alpha is not found, the value 0.1 is used as fallback.
   * @param[in]     mesh       @copybrief P3MParameters::mesh
   * @param[in]     cao        @copybrief P3MParameters::cao
   * @param[in]     r_cut_iL   @copybrief P3MParameters::r_cut_iL
   * @returns Error magnitude, real-space error, k-space error,
   *          @copybrief P3MParameters::alpha_L
   */
  virtual std::tuple<double, double, double, double>
  calculate_accuracy(Utils::Vector3i const &mesh, int cao,
                     double r_cut_iL) const = 0;

  /** @brief Veto real-space cutoffs larger than the layer correction gap. */
  virtual std::optional<std::string>
  layer_correction_veto_r_cut(double r_cut) const = 0;

  /** @brief Write tuned parameters to the P3M parameter struct. */
  void commit(Utils::Vector3i const &mesh, int cao, double r_cut_iL,
              double alpha_L);

  void tune() {
    // activate tuning mode
    get_params().tuning = true;

    auto const tuned_params = get_time();

    // deactivate tuning mode
    get_params().tuning = false;

    if (tuned_params.time == time_sentinel) {
      throw std::runtime_error(m_logger->get_name() +
                               ": failed to reach requested accuracy");
    }
    // set tuned parameters
    get_params().accuracy = tuned_params.accuracy;
    commit(tuned_params.mesh, tuned_params.cao, tuned_params.r_cut_iL,
           tuned_params.alpha_L);

    m_logger->tuning_results(tuned_params.mesh, tuned_params.cao,
                             tuned_params.r_cut_iL, tuned_params.alpha_L,
                             tuned_params.accuracy, tuned_params.time);
  }

protected:
  auto get_n_trials() { return m_n_trials; }
  void increment_n_trials() { ++m_n_trials; }
  void reset_n_trials() { m_n_trials = 0ul; }
  double get_m_time(Utils::Vector3i const &mesh, int &tuned_cao,
                    double &tuned_r_cut_iL, double &tuned_alpha_L,
                    double &tuned_accuracy);
  double get_mc_time(Utils::Vector3i const &mesh, int cao,
                     double &tuned_r_cut_iL, double &tuned_alpha_L,
                     double &tuned_accuracy);
};

#endif // P3M or DP3M
