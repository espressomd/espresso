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

#pragma once

/**
 * @file caliper_utils.hpp
 * @brief Zero-overhead Caliper guards for the inactive (no @c CALI_CONFIG)
 * case.
 *
 * When Caliper is compiled in but @c CALI_CONFIG is not set (the common
 * production case), the standard @c CALI_MARK_BEGIN/END and
 * @c CALI_CXX_MARK_FUNCTION macros still pay the full Caliper entry cost
 * (thread-blackboard update, siglock acquisition).  On a hot path with
 * many markers per step this amounts to measurable wall-clock overhead.
 *
 * This header provides guarded replacements:
 * - @c espresso_cali_active() — checked once per process; in steady state this
 *   is a single byte load + branch.
 * - @c ESPRESSO_CALI_MARK_FUNCTION — RAII guard equivalent to
 *   @c CALI_CXX_MARK_FUNCTION but zero-cost when inactive.
 * - @c ESPRESSO_CALI_MARK_BEGIN(name) / @c ESPRESSO_CALI_MARK_END(name)
 * - @c EspressoCaliLoop — RAII loop wrapper replacing the
 *   @c ESPRESSO_CALI_MARK_LOOP_BEGIN / @c ESPRESSO_CALI_MARK_LOOP_ITERATION /
 *   @c ESPRESSO_CALI_MARK_LOOP_END macro triplet.
 *
 * All macros and types are no-ops / empty when @c ESPRESSO_CALIPER is not
 * defined.
 *
 * @par Activation path limitation
 * @c espresso_cali_active() detects Caliper activation **solely** via the
 * @c CALI_CONFIG environment variable.  Activation through
 * @c CALI_SERVICES_ENABLE, Caliper config files, or programmatic
 * @c cali::ConfigManager::start() is **not** detected — hot-path regions will
 * silently produce no output for those activation paths.  @c CALI_CONFIG is
 * the supported activation path for ESPResSo profiling (used by
 * @c testsuite/python/caliper.py).  This restriction is intentional: checking
 * @c CALI_CONFIG once and caching the result eliminates the Caliper entry cost
 * (siglock + thread-blackboard update) on every hot-path marker when profiling
 * is inactive.
 */

#include <config/config.hpp>

#ifdef ESPRESSO_CALIPER

#include <caliper/cali.h>
#include <cstdlib>

/**
 * @brief Return true if Caliper is configured for this process.
 *
 * Reads @c CALI_CONFIG from the environment exactly once (on the first call)
 * and caches the result.  Subsequent calls pay only one byte load + branch.
 *
 * @c inline ensures a single shared static across all translation units
 * (C++ ODR for inline functions with static locals).
 *
 * @note Only @c CALI_CONFIG activation is detected; see the file-level
 *       documentation for the rationale and limitation.
 */
inline bool espresso_cali_active() noexcept {
  static const bool active = (std::getenv("CALI_CONFIG") != nullptr);
  return active;
}

/**
 * @brief RAII region guard: begin on construction, end on destruction.
 *
 * Conditionally calls @c cali_begin_region / @c cali_end_region only when
 * @c espresso_cali_active() is true.  Correctly handles all exit paths.
 */
struct EspressoCaliRegion {
  const char *name_;
  bool active_;

  EspressoCaliRegion(const char *name, bool active)
      : name_(name), active_(active) {
    if (active_)
      cali_begin_region(name_);
  }
  ~EspressoCaliRegion() {
    if (active_)
      cali_end_region(name_);
  }
  EspressoCaliRegion(const EspressoCaliRegion &) = delete;
  EspressoCaliRegion &operator=(const EspressoCaliRegion &) = delete;
};

/**
 * @brief RAII iteration annotation for one pass through a loop body.
 *
 * Returned by @c EspressoCaliLoop::iteration().  When active, calls
 * @c cali_begin_int on the loop's iteration attribute on construction and
 * @c cali_end_byid on destruction — identical to what @c cali::Loop::Iteration
 * does internally.  When inactive, both operations are no-ops.
 *
 * The object must be stored at the top of the loop body so that its lifetime
 * spans the loop body:
 * @code
 *   for (int step = 0; ...) {
 *     auto cali_iter = cali_loop.iteration(step);
 *     // ...
 *   }  // cali_iter destroyed here, closing the iteration region
 * @endcode
 */
struct EspressoCaliIteration {
  cali_id_t iter_attr_;
  bool active_;

  EspressoCaliIteration(cali_id_t iter_attr, int iter, bool active)
      : iter_attr_(iter_attr), active_(active) {
    if (active_)
      cali_begin_int(iter_attr_, iter);
  }
  ~EspressoCaliIteration() {
    if (active_)
      cali_end(iter_attr_);
  }
  EspressoCaliIteration(const EspressoCaliIteration &) = delete;
  EspressoCaliIteration &operator=(const EspressoCaliIteration &) = delete;
  // Moveable so it can be returned from iteration() with NRVO/move.
  EspressoCaliIteration(EspressoCaliIteration &&other) noexcept
      : iter_attr_(other.iter_attr_), active_(other.active_) {
    other.active_ = false;
  }
  EspressoCaliIteration &operator=(EspressoCaliIteration &&) = delete;
};

/**
 * @brief RAII loop wrapper replacing the @c CALI_CXX_MARK_LOOP_BEGIN /
 * @c CALI_CXX_MARK_LOOP_ITERATION / @c CALI_CXX_MARK_LOOP_END macro triplet.
 *
 * When @c espresso_cali_active() is true at construction time, creates the
 * same @c cali.loop region and @c iteration#name attribute as the raw Caliper
 * macros would.  When inactive, construction, iteration(), and destruction are
 * all no-ops with no Caliper calls.
 *
 * Usage:
 * @code
 *   EspressoCaliLoop cali_loop("Integration loop");
 *   for (int step = 0; step < n_steps; ++step) {
 *     auto cali_iter = cali_loop.iteration(step);  // RAII: closed at }
 *     // ... loop body ...
 *   }
 *   // cali_loop destructor ends the loop region
 * @endcode
 */
struct EspressoCaliLoop {
  cali::Loop *loop_;
  cali_id_t iter_attr_;

  explicit EspressoCaliLoop(const char *name)
      : loop_(nullptr), iter_attr_(CALI_INV_ID) {
    if (espresso_cali_active()) {
      loop_ = new cali::Loop(name);
      iter_attr_ = cali_make_loop_iteration_attribute(name);
    }
  }
  ~EspressoCaliLoop() {
    if (loop_)
      loop_->end();
    delete loop_;
  }
  EspressoCaliLoop(const EspressoCaliLoop &) = delete;
  EspressoCaliLoop &operator=(const EspressoCaliLoop &) = delete;

  /**
   * @brief Return an RAII iteration annotation for the current step.
   *
   * The returned @c EspressoCaliIteration must be kept alive for the duration
   * of the loop body; store it in a named local variable.
   */
  EspressoCaliIteration iteration(int iter) const {
    return EspressoCaliIteration(iter_attr_, iter, loop_ != nullptr);
  }
};

// NOLINTBEGIN(cppcoreguidelines-macro-usage)

/**
 * @brief Guarded drop-in replacement for @c CALI_CXX_MARK_FUNCTION.
 *
 * Captures @c __func__ at the call site so that the region name is the
 * enclosing function name, not the destructor name.
 */
#define ESPRESSO_CALI_MARK_FUNCTION                                            \
  EspressoCaliRegion CALI_CREATE_VAR_NAME(__espresso_cali_fn, __LINE__)(       \
      __func__, ::espresso_cali_active())

/**
 * @brief Guarded @c CALI_MARK_BEGIN — no-op when inactive.
 */
#define ESPRESSO_CALI_MARK_BEGIN(name)                                         \
  do {                                                                         \
    if (::espresso_cali_active())                                              \
      CALI_MARK_BEGIN(name);                                                   \
  } while (false)

/**
 * @brief Guarded @c CALI_MARK_END — no-op when inactive.
 */
#define ESPRESSO_CALI_MARK_END(name)                                           \
  do {                                                                         \
    if (::espresso_cali_active())                                              \
      CALI_MARK_END(name);                                                     \
  } while (false)

// NOLINTEND(cppcoreguidelines-macro-usage)

#endif // ESPRESSO_CALIPER
