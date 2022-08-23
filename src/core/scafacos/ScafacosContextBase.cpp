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

#include "config/config.hpp"

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLAR)

#include "scafacos/ScafacosContextBase.hpp"

#include <boost/algorithm/string/join.hpp>

#include <scafacos/Scafacos.hpp>

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::string> ScafacosContextBase::available_methods() {
  return ::Scafacos::Scafacos::available_methods();
}

void ScafacosContextBase::sanity_check_method(std::string const &method_name) {
  auto const all = ScafacosContextBase::available_methods();
  auto const valid_methods = std::set<std::string>(all.begin(), all.end());
  if (valid_methods.count(method_name) == 0) {
    auto const method_names = "'" + boost::algorithm::join(all, "', '") + "'";
    throw std::invalid_argument("Method '" + method_name +
                                "' is unknown or not compiled in ScaFaCoS; "
                                "currently available methods are " +
                                method_names);
  }
}

#endif // SCAFACOS or SCAFACOS_DIPOLAR
