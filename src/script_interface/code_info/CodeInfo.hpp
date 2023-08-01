/*
 * Copyright (C) 2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_CODE_INFO_CODE_INFO_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_CODE_INFO_CODE_INFO_HPP

#include "script_interface/ScriptInterface.hpp"

#include <string>
#include <vector>

namespace ScriptInterface {
namespace CodeInfo {

class CodeInfo : public ObjectHandle {
public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &parameters) override;
};

void check_features(std::vector<std::string> const &features);

} // namespace CodeInfo
} // namespace ScriptInterface

#endif
