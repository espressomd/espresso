/*
 * Copyright (C) 2019-2023 The ESPResSo project
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

#ifdef CALIPER

#include "script_interface/ScriptInterface.hpp"

#include <caliper/cali.h>

#include <string>

namespace ScriptInterface::Profiler {

class Caliper : public ObjectHandle {
public:
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "begin_section") {
      auto const label = get_value<std::string>(parameters, "label");
      CALI_MARK_BEGIN(label.c_str());
      return {};
    }
    if (method == "end_section") {
      auto const label = get_value<std::string>(parameters, "label");
      CALI_MARK_END(label.c_str());
      return {};
    }
    return {};
  }
};

} // namespace ScriptInterface::Profiler

#endif // CALIPER
