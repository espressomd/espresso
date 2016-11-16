#ifndef SCRIPT_INTERFACE_VARIANT_HPP
#define SCRIPT_INTERFACE_VARIANT_HPP

#include <boost/variant.hpp>

#include "core/Vector.hpp"
#include "utils/ObjectId.hpp"

namespace ScriptInterface {
class ScriptInterfaceBase;
using ObjectId = Utils::ObjectId<ScriptInterfaceBase>;

/**
 * @brief Possible types for parameters.
 */

typedef boost::make_recursive_variant<
    bool, int, double, std::string, std::vector<int>, std::vector<double>,
    ObjectId, std::vector<boost::recursive_variant_>>::type Variant;

enum class VariantType {
  BOOL = 0,
  INT,
  DOUBLE,
  STRING,
  INT_VECTOR,
  DOUBLE_VECTOR,
  OBJECTID,
  VECTOR
};

std::string get_type_label(Variant const &);
std::string get_type_label(VariantType);

bool is_bool(Variant const &v);
bool is_int(Variant const &v);
bool is_string(Variant const &v);
bool is_double(Variant const &v);
bool is_int_vector(Variant const &v);
bool is_double_vector(Variant const &v);
bool is_objectid(Variant const &v);
bool is_vector(Variant const &v);

void transform_vectors(Variant &v);

std::string print_variant_types(Variant const &v);

} /* namespace ScriptInterface */

#endif
