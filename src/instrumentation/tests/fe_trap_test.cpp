/*
 * Copyright (C) 2024 The ESPResSo project
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

#define BOOST_TEST_MODULE "floating-point exceptions test"
#define BOOST_TEST_DYN_LINK

#include <cfenv>

#if defined(__STDC_IEC_559__) and defined(__GLIBC__) and defined(__x86_64__)
#define ESPRESSO_FPE_IS_SUPPORTED
#define ESPRESSO_FPE_USING_GLIBC_X86_64
#elif defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FPE_IS_SUPPORTED
#define ESPRESSO_FPE_USING_APPLE_ARM_64
#endif

#if defined(ESPRESSO_FPE_IS_SUPPORTED) and not defined(__FAST_MATH__)
#include <boost/test/unit_test.hpp>

#include <instrumentation/fe_trap.hpp>

#include <stdexcept>

#if defined(ESPRESSO_FPE_USING_GLIBC_X86_64)
#include <cassert>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <initializer_list>

static volatile std::sig_atomic_t last_signal_status = 0;
static volatile std::sig_atomic_t last_signal_code = 0;
static std::jmp_buf jmp_env;

static void fpe_signal_handler(int signum, siginfo_t *sip, void *) {
  /* Signal handlers always return to the instruction that failed,
   * which leads to an infinite loop if the signal was sent by the
   * instruction itself.
   *
   * In the case of @c SIGFPE, we cannot design an abnormal mathematical
   * operation that relies on a volatile-qualified variable that can be
   * changed from within the signal handler (undefined behavior).
   * But we can long jump out of the signal handler,
   * as long as the signal wasn't sent by an operation that isn't
   * async-signal-safe, see POSIX.1-2008 Technical Corrigendum 2:
   * https://www.man7.org/linux/man-pages/man3/sigsetjmp.3.html#:~:text=Undefined%20Behavior
   * and the GNU libc relevant page:
   * https://www.gnu.org/savannah-checkouts/gnu/libc/manual/html_node/Longjmp-in-Handler.html
   * A call to @c siglongjmp is needed to preserve the signal mask, otherwise
   * the signal handler won't be able to catch any more signals, see:
   * https://pubs.opengroup.org/onlinepubs/9699919799/functions/sigaction.html#tag_16_540_07
   */
  ::last_signal_status = signum;
  ::last_signal_code = sip->si_code;
  siglongjmp(::jmp_env, 1);
}

static void set_sigaction(int signum, struct sigaction const *new_action,
                          struct sigaction *old_action) {
  [[maybe_unused]] auto const ret = sigaction(signum, new_action, old_action);
  assert(ret == 0);
}

BOOST_AUTO_TEST_CASE(trap_by_signal) {
  // declare local variables volatile to avoid compiler optimizations
  double volatile bad_denominator = 0.;
  double volatile bad_exponent = 10000.;
  double volatile value = 1.;
  // without instrumentation, abnormal operations are allowed
  value = std::exp(bad_exponent);
  value = 0. / bad_denominator;
  BOOST_REQUIRE(std::isnan(value));
  BOOST_REQUIRE_EQUAL(::last_signal_status, 0);
  BOOST_REQUIRE_EQUAL(::last_signal_code, 0);

  // catch signals sent by FPE traps
  struct sigaction old_action;
  struct sigaction fpe_action;
  fpe_action.sa_sigaction = &fpe_signal_handler;
  sigemptyset(&fpe_action.sa_mask);
  fpe_action.sa_flags = SA_SIGINFO;
  set_sigaction(SIGFPE, &fpe_action, &old_action);

  {
    auto const trap = fe_trap::make_unique_scoped();
    BOOST_REQUIRE(trap.is_unique());
    value = 0.;
    while (sigsetjmp(::jmp_env, 1) == 0) {
      // LCOV_EXCL_START
      value = 2.;
      value = 0. / bad_denominator;
      // LCOV_EXCL_STOP
    }
    BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
    BOOST_CHECK_EQUAL(::last_signal_code, FPE_FLTINV);
    BOOST_REQUIRE(not std::isnan(value));
    BOOST_REQUIRE_EQUAL(value, 2.);
    ::last_signal_status = 0;
    ::last_signal_code = 0;
  }
  {
    auto const trap = fe_trap::make_shared_scoped();
    BOOST_REQUIRE(not trap.is_unique());
    value = 0.;
    while (sigsetjmp(::jmp_env, 1) == 0) {
      // LCOV_EXCL_START
      value = 2.;
      value = std::exp(bad_exponent);
      // LCOV_EXCL_STOP
    }
    BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
    BOOST_CHECK_EQUAL(::last_signal_code, FPE_FLTOVF);
    BOOST_REQUIRE(not std::isnan(value));
    BOOST_REQUIRE_EQUAL(value, 2.);
    ::last_signal_status = 0;
    ::last_signal_code = 0;
  }
  {
    auto const trap1 = fe_trap::make_shared_scoped(FE_UNDERFLOW);
    {
      auto const trap2 = fe_trap::make_shared_scoped(FE_UNDERFLOW);
      BOOST_REQUIRE_EQUAL(trap1.get_flags(), trap2.get_flags());
      value = 0.;
      while (sigsetjmp(::jmp_env, 1) == 0) {
        // LCOV_EXCL_START
        value = 2.;
        value = std::exp(-bad_exponent);
        // LCOV_EXCL_STOP
      }
      BOOST_CHECK_EQUAL(::last_signal_status, SIGFPE);
      BOOST_CHECK_EQUAL(::last_signal_code, FPE_FLTUND);
      BOOST_REQUIRE(not std::isnan(value));
      BOOST_REQUIRE_EQUAL(value, 2.);
      ::last_signal_status = 0;
      ::last_signal_code = 0;
    }
  }

  // reset default signal handler
  set_sigaction(SIGFPE, &old_action, nullptr);
  // without instrumentation, abnormal operations are allowed
  value = 1. / bad_denominator;
  value = std::exp(bad_exponent);
}
#endif // defined(ESPRESSO_FPE_USING_GLIBC_X86_64)

BOOST_AUTO_TEST_CASE(exceptions) {
  {
    auto const trap = fe_trap::make_unique_scoped();
    BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
    BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(), std::runtime_error);
  }
  {
    auto const trap = fe_trap::make_shared_scoped();
    // cannot instantiate a unique object when a shared object already exists
    BOOST_REQUIRE_THROW(fe_trap::make_unique_scoped(), std::runtime_error);
    // cannot instantiate a shared object with mismatched flags (singleton)
    for (int other_excepts : {FE_INEXACT, FE_ALL_EXCEPT}) {
      if (trap.get_flags() != other_excepts) {
        BOOST_REQUIRE_THROW(fe_trap::make_shared_scoped(FE_ALL_EXCEPT),
                            std::invalid_argument);
      }
    }
  }
}
#else
int main() {}
#endif // defined(ESPRESSO_FPE_IS_SUPPORTED) and not defined(__FAST_MATH__)
