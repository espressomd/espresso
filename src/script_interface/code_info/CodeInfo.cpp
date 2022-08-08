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

#include "CodeInfo.hpp"

#include "config-features.hpp"
#include "script_interface/scafacos/scafacos.hpp"
#include "version.hpp"

#include <string>
#include <vector>

namespace ScriptInterface {
namespace CodeInfo {

static Variant get_feature_list(char const *const ptr[], unsigned int len) {
  return make_vector_of_variants(std::vector<std::string>{ptr, ptr + len});
}

Variant CodeInfo::do_call_method(std::string const &name,
                                 VariantMap const &parameters) {
  if (name == "features") {
    return get_feature_list(FEATURES, NUM_FEATURES);
  }
  if (name == "all_features") {
    return get_feature_list(FEATURES_ALL, NUM_FEATURES_ALL);
  }
  if (name == "build_type") {
    return std::string(ESPRESSO_BUILD_TYPE);
  }
  if (name == "scafacos_methods") {
#ifdef SCAFACOS
    return make_vector_of_variants(Scafacos::available_methods());
#else  // SCAFACOS
    return make_vector_of_variants(std::vector<std::string>(0));
#endif // SCAFACOS
  }
  return {};
}

} // namespace CodeInfo
} // namespace ScriptInterface
