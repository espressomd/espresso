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
#ifndef PROFILER_PROFILER_HPP
#define PROFILER_PROFILER_HPP

#include <string>

#ifdef HAVE_CALIPER
#include <caliper/cali.h>

#define ESPRESSO_PROFILER_CXX_MARK_FUNCTION CALI_CXX_MARK_FUNCTION
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN CALI_CXX_MARK_LOOP_BEGIN
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_END CALI_CXX_MARK_LOOP_END
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION CALI_CXX_MARK_LOOP_ITERATION
#define ESPRESSO_PROFILER_MARK_FUNCTION_BEGIN CALI_MARK_FUNCTION_BEGIN
#define ESPRESSO_PROFILER_MARK_FUNCTION_END CALI_MARK_FUNCTION_END
#define ESPRESSO_PROFILER_MARK_LOOP_BEGIN CALI_MARK_LOOP_BEGIN
#define ESPRESSO_PROFILER_MARK_LOOP_END CALI_MARK_LOOP_END
#define ESPRESSO_PROFILER_MARK_ITERATION_BEGIN CALI_MARK_ITERATION_BEGIN
#define ESPRESSO_PROFILER_MARK_ITERATION_END CALI_MARK_ITERATION_END
#define ESPRESSO_PROFILER_WRAP_STATEMENT CALI_WRAP_STATEMENT
#define ESPRESSO_PROFILER_MARK_BEGIN CALI_MARK_BEGIN
#define ESPRESSO_PROFILER_MARK_END CALI_MARK_END
#else
#define ESPRESSO_PROFILER_CXX_MARK_FUNCTION
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN(A, B)
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_END(A)
#define ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION(A, B)
#define ESPRESSO_PROFILER_MARK_FUNCTION_BEGIN
#define ESPRESSO_PROFILER_MARK_FUNCTION_END
#define ESPRESSO_PROFILER_MARK_LOOP_BEGIN(A, B)
#define ESPRESSO_PROFILER_MARK_LOOP_END(A)
#define ESPRESSO_PROFILER_MARK_ITERATION_BEGIN(A, B)
#define ESPRESSO_PROFILER_MARK_ITERATION_END(A)
#define ESPRESSO_PROFILER_WRAP_STATEMENT(A, B)
#define ESPRESSO_PROFILER_MARK_BEGIN(A)
#define ESPRESSO_PROFILER_MARK_END(A)
#endif

namespace Profiler {
/**
 * @brief Start named section.
 *
 * @param name Identifier of the section.
 */
inline void begin_section(const std::string &name) {
  ESPRESSO_PROFILER_MARK_BEGIN(name.c_str());
}

/**
 * @brief End named section.
 *
 * @param name Identifier of the section.
 */
inline void end_section(const std::string &name) {
  ESPRESSO_PROFILER_MARK_BEGIN(name.c_str());
}
} // namespace Profiler
#endif
