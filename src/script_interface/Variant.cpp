/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Variant.hpp"

namespace ScriptInterface {
static const char *VariantLabels[] = {"NONE",          "BOOL",     "INT",
                                      "DOUBLE",        "STRING",   "INT_VECTOR",
                                      "DOUBLE_VECTOR", "OBJECTID", "VECTOR"};

std::string get_type_label(Variant const &v) {
  return std::string(VariantLabels[v.which()]);
}

std::string get_type_label(VariantType t) {
  return std::string(VariantLabels[static_cast<int>(t)]);
}

namespace {
template <typename T>
std::vector<T> to_vector(std::vector<Variant> const &variant_vector) {
  std::vector<T> ret;
  ret.reserve(variant_vector.size());

  for (auto const &it : variant_vector) {
    ret.emplace_back(boost::get<T>(it));
  }

  return ret;
}
} /* namespace */

void transform_vectors(Variant &v) {
  if (is_type<std::vector<Variant>>(v)) {
    auto &variant_vector = boost::get<std::vector<Variant>>(v);

    /* only int, transform to vector<int> */
    if (std::all_of(variant_vector.begin(), variant_vector.end(), is_type<int>)) {
      v = to_vector<int>(variant_vector);
      return;
    }

    /* only double, transform to vector<double> */
    if (std::all_of(variant_vector.begin(), variant_vector.end(), is_type<double>)) {
      v = to_vector<double>(variant_vector);
      return;
    }

    /* v is a mixed vector, recurse into the elements. */
    for (auto &it : variant_vector) {
      transform_vectors(it);
    }
  }
}

std::string print_variant_types(Variant const &v) {
  if (is_type<std::vector<Variant>>(v)) {
    auto const &variant_vector = boost::get<std::vector<Variant>>(v);
    std::string ret{"{"};

    for (auto const &it : variant_vector) {
      ret += print_variant_types(it);
      ret += ", ";
    }
    ret += "}";

    return ret;
  } else {
    return get_type_label(v);
  }
}

} /* namespace ScriptInterface */
