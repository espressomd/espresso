/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#ifdef ESPRESSO_FPE

#include <memory>
#include <mutex>
#include <optional>
#include <utility>

/**
 * @brief Floating-point exception trap.
 *
 * Thread-safe RAII-style mechanism to trap floating-point exceptions
 * (on platforms that support floating-point environment management) for
 * the duration of a scoped block. Exception traps are set when the object
 * is created; when getting out-of-scope, either normally or during stack
 * unwinding, the exception traps are automatically reset.
 *
 * Please note "exception" and "exception handling" have a specific meaning
 * in this context and are completely unrelated to C++ exceptions.
 * For more details, see annex F IEC 60559 "floating-point arithmetic"
 * in ISO/EIC 9899 @cite ISO-EIC-9899-1999 and chapter 7
 * "Exceptions and default exception handling" in
 * ISO/IEC 60559:2020(E) @cite ISO-EIC-60559-2020.
 *
 * The exception handling behavior is implementation-defined. For example,
 * GNU libc sends the @c SIGFPE signal on x86 architectures; it can be caught
 * by a signal handler that leverages stack environments and long jumps via
 * [`sigsetjmp()`](https://www.man7.org/linux/man-pages/man3/sigsetjmp.3.html).
 * On Armv8, trapping is controlled by FPCR flags; for more details,
 * see section C5.2.8 "FPCR, Floating-point Control Register" in the Armv8
 * manual @cite ARM-DDI-0487-2024. AppleClang sends the @c SIGILL signal
 * on Apple Silicon architectures.
 *
 * A modified singleton pattern is leveraged to guarantee only one trap is
 * active at any time; once expired, a new trap can be instantiated.
 * The @ref make_unique_scoped function returns a wrapper object whose
 * lifetime determines the trap duration. To set a trap in a recursive
 * function, use @ref make_shared_scoped instead.
 *
 * Usage:
 * @code{.cpp}
 * #include <instrumentation/fe_trap.hpp>
 * #include <cmath>
 * int main() {
 *   auto volatile zero = 0.;
 *   auto value = 1. / zero;  // generate NaN
 *   {
 *     auto trap = fe_trap::make_unique_scoped();
 *     value = 1. / zero;  // execution flow should be interrupted here
 *   }
 *   value = 1. / zero;  // generate NaN
 *   return std::isnan(value) ? 0 : 1;
 * }
 * @endcode
 * Build the code without fast-math and without any optimization
 * (optimizations always assume divisions by zero cannot happen):
 * @code{.sh}
 * g++ -std=c++20 -O0 -fno-fast-math -I../src/instrumentation/include \
 *     -I../src/config/include -Isrc/config/include \
 *     main.cpp ../src/instrumentation/src/fe_trap.cpp
 * ./a.out
 * }
 * @endcode
 */
class fe_trap {
  struct global_state_params {
    std::weak_ptr<fe_trap> observer;
    std::mutex mutex;
  };
  static global_state_params global_state;

  struct scoped_instance {
    explicit scoped_instance(std::shared_ptr<fe_trap> ptr)
        : m_resource{std::move(ptr)} {}
    scoped_instance(scoped_instance const &) = delete;
    scoped_instance(scoped_instance &&) noexcept = default;
    scoped_instance &operator=(scoped_instance const &) = delete;
    scoped_instance &operator=(scoped_instance &&) noexcept = default;
    bool is_unique() const { return m_resource->is_unique(); }
    int get_flags() const { return m_resource->get_flags(); }

  private:
    std::shared_ptr<fe_trap> m_resource;
  };

  struct deleter {
    void operator()(fe_trap *ptr) { delete ptr; }
  };
  friend deleter;

  int m_flags;
  bool m_unique;

  fe_trap(std::optional<int> excepts, bool unique);
  ~fe_trap();

  static int parse_excepts(std::optional<int> excepts);

public:
  fe_trap(fe_trap const &) = delete;
  fe_trap(fe_trap &&) noexcept = delete;
  fe_trap &operator=(fe_trap const &) = delete;
  fe_trap &operator=(fe_trap &&) noexcept = delete;
  /** @brief Get floating-point exception flags. */
  int get_flags() const { return m_flags; }
  /** @brief Check if this handle is a unique handle. */
  bool is_unique() const { return m_unique; }

  /**
   * @brief Generate a unique trap with the lifetime of the current scope.
   * @param excepts Combination of floating-point exception flags.
   */
  static scoped_instance
  make_unique_scoped(std::optional<int> excepts = std::nullopt);
  /**
   * @brief Generate a shared trap with the lifetime of the current scope.
   * Subsequent calls to this function will yield the same trap handle,
   * as long as they have the same parameter @c excepts.
   * @param excepts Combination of floating-point exception flags.
   */
  static scoped_instance
  make_shared_scoped(std::optional<int> excepts = std::nullopt);
};

#endif // ESPRESSO_FPE
