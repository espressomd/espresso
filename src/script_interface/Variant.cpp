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

bool is_none(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::NONE);
}

bool is_bool(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::BOOL);
}

bool is_int(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::INT);
}

bool is_string(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::STRING);
}

bool is_double(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::DOUBLE);
}

bool is_int_vector(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::INT_VECTOR);
}

bool is_double_vector(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::DOUBLE_VECTOR);
}

bool is_objectid(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::OBJECTID);
}

bool is_vector(Variant const &v) {
  return v.which() == static_cast<int>(VariantType::VECTOR);
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
  if (is_vector(v)) {
    auto &variant_vector = boost::get<std::vector<Variant>>(v);

    /* only int, tranform to vector<int> */
    if (std::all_of(variant_vector.begin(), variant_vector.end(), is_int)) {
      v = to_vector<int>(variant_vector);
      return;
    }

    /* only double, tranform to vector<int> */
    if (std::all_of(variant_vector.begin(), variant_vector.end(), is_double)) {
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
  if (is_vector(v)) {
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
