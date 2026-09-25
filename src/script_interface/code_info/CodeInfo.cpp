/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include "config/config-features.hpp"
#include "config/config-features.impl.hpp"
#include "config/version.hpp"
#include "script_interface/scafacos/scafacos.hpp"

#include <boost/algorithm/string/join.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace ScriptInterface {
namespace CodeInfo {

static auto get_feature_vector(char const *const ptr[], std::size_t len) {
  return std::vector<std::string>{ptr, ptr + len};
}

static auto get_feature_set(char const *const ptr[], std::size_t len) {
  return std::unordered_set<std::string>(ptr, ptr + len);
}

Variant CodeInfo::do_call_method(std::string const &name, VariantMap const &) {
  if (name == "features") {
    return make_vector_of_variants(get_feature_vector(FEATURES, NUM_FEATURES));
  }
  if (name == "all_features") {
    return make_vector_of_variants(
        get_feature_vector(FEATURES_ALL, NUM_FEATURES_ALL));
  }
  if (name == "build_type") {
    return std::string(ESPRESSO_BUILD_TYPE);
  }
  if (name == "scafacos_methods") {
#ifdef ESPRESSO_SCAFACOS
    return make_vector_of_variants(Scafacos::available_methods());
#else  // ESPRESSO_SCAFACOS
    return make_vector_of_variants(std::vector<std::string>(0));
#endif // ESPRESSO_SCAFACOS
  }
  if (name == "has_fast_math") {
#if defined(__FAST_MATH__)
    return true;
#else
    return false;
#endif
  }
  if (name == "toolchain") {
#define ESPRESSO_REGISTER_METADATA(name)                                       \
  flat_map[#name] = serializer(ESPRESSO_##name);
    VariantMap flat_map;
    auto const serializer = [](std::string const value) {
      Variant result = None{};
      if (not value.empty()) {
        result = value;
      }
      return result;
    };
    ESPRESSO_REGISTER_METADATA(CMAKE_C_COMPILER_ID)
    ESPRESSO_REGISTER_METADATA(CMAKE_C_COMPILER_VERSION)
    ESPRESSO_REGISTER_METADATA(CMAKE_CXX_COMPILER_ID)
    ESPRESSO_REGISTER_METADATA(CMAKE_CXX_COMPILER_VERSION)
    ESPRESSO_REGISTER_METADATA(CMAKE_CUDA_COMPILER_ID)
    ESPRESSO_REGISTER_METADATA(CMAKE_CUDA_COMPILER_VERSION)
    ESPRESSO_REGISTER_METADATA(CMAKE_CUDA_HOST_COMPILER_ID)
    ESPRESSO_REGISTER_METADATA(CMAKE_CUDA_HOST_COMPILER_VERSION)
    ESPRESSO_REGISTER_METADATA(OpenMP_VERSION)
    ESPRESSO_REGISTER_METADATA(OpenMP_C_VERSION)
    ESPRESSO_REGISTER_METADATA(OpenMP_CXX_VERSION)
    ESPRESSO_REGISTER_METADATA(OpenMP_CUDA_VERSION)
    ESPRESSO_REGISTER_METADATA(ESPRESSO_MPIEXEC_VENDOR)
    ESPRESSO_REGISTER_METADATA(ESPRESSO_MPIEXEC_VERSION)
    return flat_map;
  }
  return {};
}

std::string check_features_msg(std::vector<std::string> const &features) {
  auto const allowed = get_feature_set(FEATURES_ALL, NUM_FEATURES_ALL);
  auto const compiled_features = get_feature_set(FEATURES, NUM_FEATURES);
  std::vector<std::string> missing_features{};
  for (auto const &feature : features) {
    if (not allowed.contains(feature)) {
      return "Unknown feature '" + feature + "'";
    }
    if (not compiled_features.contains(feature)) {
      missing_features.emplace_back(feature);
    }
  }
  if (not missing_features.empty()) {
    return "Missing features " + boost::algorithm::join(missing_features, ", ");
  }
  return "";
}

void check_features(std::vector<std::string> const &features) {
  auto const error_msg = check_features_msg(features);
  if (not error_msg.empty()) {
    throw std::runtime_error(error_msg);
  }
}

} // namespace CodeInfo
} // namespace ScriptInterface
