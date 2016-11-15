#include "Variant.hpp"

namespace ScriptInterface {
const char *VariantLabels[] = {
    "BOOL",          "INT",      "DOUBLE",  "STRING",   "INT_VECTOR",
    "DOUBLE_VECTOR", "VECTOR2D", "VECTOR3", "OBJECTID", "VECTOR"};

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
} /* namespace ScriptInterface */
