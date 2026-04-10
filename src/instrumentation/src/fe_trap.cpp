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

#include <config/config.hpp>

#ifdef ESPRESSO_FPE

#include <instrumentation/fe_trap.hpp>

#include <cassert>
#include <cfenv>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>

#if defined(__STDC_IEC_559__) and defined(__GLIBC__) and defined(__x86_64__)
#define ESPRESSO_FPE_USING_GLIBC_X86_64
#elif defined(__arm64__) and defined(__APPLE__)
#define ESPRESSO_FPE_USING_APPLE_ARM_64
#endif

fe_trap::global_state_params fe_trap::global_state{{}, {}};

void fe_trap::activate() {
  if (m_active) {
    return;
  }
#if defined(ESPRESSO_FPE_USING_GLIBC_X86_64)
  {
    [[maybe_unused]] auto const status = feenableexcept(m_flags);
    // note: status should be 0 since we use the singleton pattern
    assert(status == 0);
  }
#elif defined(ESPRESSO_FPE_USING_APPLE_ARM_64)
  {
    using fpcr_t = decltype(std::fenv_t::__fpcr);
    std::fenv_t env;
    {
      [[maybe_unused]] auto const status = std::fegetenv(&env);
      assert(status == 0u);
    }
    env.__fpcr |= static_cast<fpcr_t>(m_flags);
    {
      [[maybe_unused]] auto const status = std::fesetenv(&env);
      assert(status == 0u);
    }
  }
#else
#error "FE not supported"
#endif
  m_active = true;
}

void fe_trap::deactivate() {
  if (not m_active) {
    return;
  }
#if defined(ESPRESSO_FPE_USING_GLIBC_X86_64)
  {
    [[maybe_unused]] auto const status = fedisableexcept(m_flags);
    // note: status can become 0 in a signal handler that calls @c siglongjmp()
    assert(status == 0 or status == m_flags);
  }
#elif defined(ESPRESSO_FPE_USING_APPLE_ARM_64)
  {
    using fpcr_t = decltype(std::fenv_t::__fpcr);
    std::fenv_t env;
    {
      [[maybe_unused]] auto const status = std::fegetenv(&env);
      assert(status == 0u);
    }
    assert((env.__fpcr & static_cast<fpcr_t>(m_flags)) ==
           static_cast<fpcr_t>(m_flags));
    env.__fpcr &= static_cast<fpcr_t>(~m_flags);
    {
      [[maybe_unused]] auto const status = std::fesetenv(&env);
      assert(status == 0u);
    }
  }
#else
#error "FE not supported"
#endif
  m_active = false;
}

int fe_trap::parse_excepts(std::optional<int> excepts) {
  auto const fallback = FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW;
  int retval = excepts ? *excepts : fallback;
#if defined(ESPRESSO_FPE_USING_APPLE_ARM_64)
  retval <<= 8;
#endif
  return retval;
}

fe_trap::scoped_instance
fe_trap::make_unique_scoped(std::optional<int> excepts) {
  std::lock_guard<std::mutex> lock(fe_trap::global_state.mutex);
  if (fe_trap::global_state.observer.lock()) {
    throw std::runtime_error("Cannot create more than 1 instance of fe_trap");
  }
  auto raw_ptr = new fe_trap(excepts, true);
  auto watched = std::shared_ptr<fe_trap>(raw_ptr, deleter{});
  fe_trap::global_state.observer = watched;
  return fe_trap::scoped_instance(watched);
}

fe_trap::scoped_instance
fe_trap::make_shared_scoped(std::optional<int> excepts) {
  std::lock_guard<std::mutex> lock(fe_trap::global_state.mutex);
  if (auto watched = fe_trap::global_state.observer.lock()) {
    if (watched->is_unique()) {
      throw std::runtime_error("Cannot create more than 1 instance of fe_trap");
    }
    if (watched->get_flags() != parse_excepts(excepts)) {
      throw std::invalid_argument(
          "Cannot mix different exceptions with fe_trap");
    }
    return fe_trap::scoped_instance(watched);
  }
  auto raw_ptr = new fe_trap(excepts, false);
  auto watched = std::shared_ptr<fe_trap>(raw_ptr, deleter{});
  fe_trap::global_state.observer = watched;
  return fe_trap::scoped_instance(watched);
}

std::shared_ptr<fe_trap::scoped_pause> fe_trap::make_shared_pause_scoped() {
  std::lock_guard<std::mutex> lock(fe_trap::global_state.mutex);
  if (auto watched = fe_trap::global_state.observer.lock()) {
    if (auto pause = watched->m_pause.lock()) {
      return pause;
    }
    auto pause = std::make_shared<scoped_pause>(watched);
    watched->m_pause = pause;
    return pause;
  }
  return {};
}

#endif // ESPRESSO_FPE
