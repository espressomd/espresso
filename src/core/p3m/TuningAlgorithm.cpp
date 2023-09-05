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

#include "config/config.hpp"

#if defined(P3M) || defined(DP3M)

#include "p3m/TuningAlgorithm.hpp"
#include "p3m/common.hpp"

#include "tuning.hpp"

#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "communication.hpp"
#include "integrate.hpp"
#include "system/System.hpp"

#include <boost/range/algorithm/min_element.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>

/** @name Error codes for tuning. */
/**@{*/
/** charge assignment order too large for mesh size */
static auto constexpr P3M_TUNE_CAO_TOO_LARGE = 1.;
/** conflict with ELC gap size */
static auto constexpr P3M_TUNE_ELC_GAP_SIZE = 2.;
/** could not achieve target accuracy */
static auto constexpr P3M_TUNE_ACCURACY_TOO_LARGE = 3.;
/**@}*/

/** @brief Precision threshold for a non-zero real-space cutoff. */
static auto constexpr P3M_RCUT_PREC = 1e-3;

void TuningAlgorithm::determine_r_cut_limits() {
  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const r_cut_iL = get_params().r_cut_iL;
  if (r_cut_iL == 0.) {
    auto const min_box_l = *boost::min_element(box_geo.length());
    auto const min_local_box_l = *boost::min_element(local_geo.length());
    m_r_cut_iL_min = 0.;
    m_r_cut_iL_max = std::min(min_local_box_l, min_box_l / 2.) - skin;
    m_r_cut_iL_min *= box_geo.length_inv()[0];
    m_r_cut_iL_max *= box_geo.length_inv()[0];
  } else {
    m_r_cut_iL_min = m_r_cut_iL_max = r_cut_iL;
    m_logger->report_fixed_r_cut_iL(r_cut_iL);
  }
}

void TuningAlgorithm::determine_cao_limits(int initial_cao) {
  assert(initial_cao >= 1 and initial_cao <= 7);
  auto const cao = get_params().cao;
  if (cao == -1) {
    cao_min = 1;
    cao_max = 7;
    cao_best = initial_cao;
  } else {
    cao_min = cao_max = cao_best = cao;
    m_logger->report_fixed_cao(cao);
  }
}

void TuningAlgorithm::commit(Utils::Vector3i const &mesh, int cao,
                             double r_cut_iL, double alpha_L) {
  auto const &box_geo = *System::get_system().box_geo;
  auto &p3m_params = get_params();
  p3m_params.r_cut = r_cut_iL * box_geo.length()[0];
  p3m_params.r_cut_iL = r_cut_iL;
  p3m_params.cao = cao;
  p3m_params.alpha_L = alpha_L;
  p3m_params.alpha = alpha_L * box_geo.length_inv()[0];
  p3m_params.mesh = mesh;
}

/**
 * @brief Get the optimal alpha and the corresponding computation time
 * for a fixed @p mesh and @p cao.
 *
 * The @p tuned_r_cut_iL is determined via a simple bisection.
 *
 * @param[in]     mesh            @copybrief P3MParameters::mesh
 * @param[in]     cao             @copybrief P3MParameters::cao
 * @param[in,out] tuned_r_cut_iL  @copybrief P3MParameters::r_cut_iL
 * @param[in,out] tuned_alpha_L   @copybrief P3MParameters::alpha_L
 * @param[in,out] tuned_accuracy  @copybrief P3MParameters::accuracy
 *
 * @returns The integration time in case of success, otherwise
 *          -@ref P3M_TUNE_ACCURACY_TOO_LARGE,
 *          -@ref P3M_TUNE_CAO_TOO_LARGE, or -@ref P3M_TUNE_ELC_GAP_SIZE
 */
double TuningAlgorithm::get_mc_time(Utils::Vector3i const &mesh, int cao,
                                    double &tuned_r_cut_iL,
                                    double &tuned_alpha_L,
                                    double &tuned_accuracy) {
  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const target_accuracy = get_params().accuracy;
  double rs_err, ks_err;
  double r_cut_iL_min = m_r_cut_iL_min;
  double r_cut_iL_max = m_r_cut_iL_max;

  /* initial checks. */
  auto const k_cut_per_dir = (static_cast<double>(cao) / 2.) *
                             Utils::hadamard_division(box_geo.length(), mesh);
  auto const k_cut = *boost::min_element(k_cut_per_dir);
  auto const min_box_l = *boost::min_element(box_geo.length());
  auto const min_local_box_l = *boost::min_element(local_geo.length());
  auto const k_cut_max = std::min(min_box_l, min_local_box_l) - skin;

  if (cao >= *boost::min_element(mesh) or k_cut >= k_cut_max) {
    m_logger->log_cao_too_large(mesh[0], cao);
    return -P3M_TUNE_CAO_TOO_LARGE;
  }

  std::tie(tuned_accuracy, rs_err, ks_err, tuned_alpha_L) =
      calculate_accuracy(mesh, cao, r_cut_iL_max);

  /* Either low and high boundary are equal (for fixed cut), or the low border
     is initially 0 and therefore
     has infinite error estimate, as required. Therefore if the high boundary
     fails, there is no possible r_cut */
  if (tuned_accuracy > target_accuracy) {
    m_logger->log_skip("accuracy not achieved", mesh[0], cao, r_cut_iL_max,
                       tuned_alpha_L, tuned_accuracy, rs_err, ks_err);
    return -P3M_TUNE_ACCURACY_TOO_LARGE;
  }

  double r_cut_iL, accuracy;
  for (;;) {
    r_cut_iL = 0.5 * (r_cut_iL_min + r_cut_iL_max);

    if (r_cut_iL_max - r_cut_iL_min < P3M_RCUT_PREC)
      break;

    /* bisection */
    std::tie(accuracy, rs_err, ks_err, tuned_alpha_L) =
        calculate_accuracy(mesh, cao, r_cut_iL);
    if (accuracy > target_accuracy)
      r_cut_iL_min = r_cut_iL;
    else
      r_cut_iL_max = r_cut_iL;
  }

  /* final result is always the upper interval boundary, since only there
   * we know that the desired minimal accuracy is obtained */
  tuned_r_cut_iL = r_cut_iL = r_cut_iL_max;

  /* if we are running P3M+ELC, check that r_cut is compatible */
  auto const r_cut = r_cut_iL * box_geo.length()[0];
  auto const veto = layer_correction_veto_r_cut(r_cut);
  if (veto) {
    m_logger->log_skip(*veto, mesh[0], cao, r_cut_iL, tuned_alpha_L,
                       tuned_accuracy, rs_err, ks_err);
    return -P3M_TUNE_ELC_GAP_SIZE;
  }

  commit(mesh, cao, r_cut_iL, tuned_alpha_L);
  on_solver_change();
  auto const int_time = benchmark_integration_step(m_timings);

  std::tie(tuned_accuracy, rs_err, ks_err, tuned_alpha_L) =
      calculate_accuracy(mesh, cao, r_cut_iL);

  m_logger->log_success(int_time, mesh[0], cao, r_cut_iL, tuned_alpha_L,
                        tuned_accuracy, rs_err, ks_err);
  increment_n_trials();
  return int_time;
}

/**
 * @brief Get the optimal alpha and the corresponding computation time
 * for a fixed @p mesh.
 *
 * @p _cao should contain an initial guess, which is then adapted by stepping
 * up and down.
 *
 * @param[in]     mesh            @copybrief P3MParameters::mesh
 * @param[in,out] tuned_cao       initial guess for the
 *                                @copybrief P3MParameters::cao
 * @param[out]    tuned_r_cut_iL  @copybrief P3MParameters::r_cut_iL
 * @param[out]    tuned_alpha_L   @copybrief P3MParameters::alpha_L
 * @param[out]    tuned_accuracy  @copybrief P3MParameters::accuracy
 *
 * @returns The integration time in case of success, otherwise
 *          -@ref P3M_TUNE_CAO_TOO_LARGE
 */
double TuningAlgorithm::get_m_time(Utils::Vector3i const &mesh, int &tuned_cao,
                                   double &tuned_r_cut_iL,
                                   double &tuned_alpha_L,
                                   double &tuned_accuracy) {
  double best_time = -1., tmp_r_cut_iL = 0., tmp_alpha_L = 0.,
         tmp_accuracy = 0.;
  /* in which direction improvement is possible. Initially, we don't know it
   * yet. */
  int final_dir = 0;
  int cao = tuned_cao;

  /* the initial step sets a timing mark. If there is no valid r_cut, we can
   * only try to increase cao to increase the obtainable precision of the far
   * formula. */
  double tmp_time;
  do {
    tmp_time = get_mc_time(mesh, cao, tmp_r_cut_iL, tmp_alpha_L, tmp_accuracy);
    /* cao is too large for this grid, but still the accuracy cannot be
     * achieved, give up */
    if (tmp_time == -P3M_TUNE_CAO_TOO_LARGE) {
      return tmp_time;
    }
    /* we have a valid time, start optimising from there */
    if (tmp_time >= 0.) {
      best_time = tmp_time;
      tuned_r_cut_iL = tmp_r_cut_iL;
      tuned_alpha_L = tmp_alpha_L;
      tuned_accuracy = tmp_accuracy;
      tuned_cao = cao;
      break;
    }
    /* the required accuracy could not be obtained, try higher caos */
    cao++;
    final_dir = 1;
  } while (cao <= cao_max);
  /* with this mesh, the required accuracy cannot be obtained. */
  if (cao > cao_max)
    return -P3M_TUNE_CAO_TOO_LARGE;

  /* at the boundaries, only the opposite direction can be used for
   * optimisation
   */
  if (cao == cao_min)
    final_dir = 1;
  else if (cao == cao_max)
    final_dir = -1;

  if (final_dir == 0) {
    /* check in which direction we can optimise. Both directions are possible */
    double dir_times[3];
    for (final_dir = -1; final_dir <= 1; final_dir += 2) {
      dir_times[final_dir + 1] = tmp_time = get_mc_time(
          mesh, cao + final_dir, tmp_r_cut_iL, tmp_alpha_L, tmp_accuracy);
      /* in this direction, we cannot optimise, since we get into precision
       * trouble */
      if (tmp_time < 0.)
        continue;

      if (tmp_time < best_time) {
        best_time = tmp_time;
        tuned_r_cut_iL = tmp_r_cut_iL;
        tuned_alpha_L = tmp_alpha_L;
        tuned_accuracy = tmp_accuracy;
        tuned_cao = cao + final_dir;
      }
    }
    /* choose the direction which was optimal, if any of the two */
    if (dir_times[0] == best_time) {
      final_dir = -1;
    } else if (dir_times[2] == best_time) {
      final_dir = 1;
    } else {
      /* no improvement in either direction, however if one is only marginally
       * worse, we can still try; down is possible and not much worse, while
       * up is either illegal or even worse */
      if ((dir_times[0] >= 0 && dir_times[0] < best_time + time_granularity) &&
          (dir_times[2] < 0 || dir_times[2] > dir_times[0]))
        final_dir = -1;
      /* same for up */
      else if ((dir_times[2] >= 0 &&
                dir_times[2] < best_time + time_granularity) &&
               (dir_times[0] < 0 || dir_times[0] > dir_times[2]))
        final_dir = 1;
      else {
        /* really no chance for optimisation */
        return best_time;
      }
    }
    /* we already checked the initial cao and its neighbor */
    cao += 2 * final_dir;
  } else {
    /* here some constraint is active, and we only checked the initial cao
     * itself */
    cao += final_dir;
  }

  /* move cao into the optimisation direction until we do not gain anymore. */
  for (; cao >= cao_min && cao <= cao_max; cao += final_dir) {
    tmp_time = get_mc_time(mesh, cao, tmp_r_cut_iL, tmp_alpha_L, tmp_accuracy);
    /* if we cannot meet the precision anymore, give up */
    if (tmp_time < 0.)
      break;

    if (tmp_time < best_time) {
      best_time = tmp_time;
      tuned_r_cut_iL = tmp_r_cut_iL;
      tuned_alpha_L = tmp_alpha_L;
      tuned_accuracy = tmp_accuracy;
      tuned_cao = cao;
    } else if (tmp_time > best_time + time_granularity) {
      /* no hope of further optimisation */
      break;
    }
  }
  return best_time;
}

#endif // P3M or DP3M
