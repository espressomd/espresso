#ifndef SCRIPT_INTERFACE_VARIANT_HPP
#define SCRIPT_INTERFACE_VARIANT_HPP

#include <boost/variant.hpp>

#include "core/Vector.hpp"
#include "utils/AutoObjectId.hpp"

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

namespace detail {
/**
 * @brief Implementation of @f infer_type.
 *
 * Helper struct is needed because particle specialization
 * of functions if not allowed.
 */
template <typename T> struct infer_type_helper {};

template <> struct infer_type_helper<bool> {
  static constexpr VariantType value{VariantType::BOOL};
};

template <> struct infer_type_helper<std::string> {
  static constexpr VariantType value{VariantType::STRING};
};

template <> struct infer_type_helper<int> {
  static constexpr VariantType value{VariantType::INT};
};

template <> struct infer_type_helper<double> {
  static constexpr VariantType value{VariantType::DOUBLE};
};

template <> struct infer_type_helper<std::vector<int>> {
  static constexpr VariantType value{VariantType::INT_VECTOR};
};

template <> struct infer_type_helper<std::vector<double>> {
  static constexpr VariantType value{VariantType::DOUBLE_VECTOR};
};

template <size_t N> struct infer_type_helper<Vector<N, double>> {
  static constexpr VariantType value{VariantType::DOUBLE_VECTOR};
};

template <size_t N> struct infer_type_helper<Vector<N, int>> {
  static constexpr VariantType value{VariantType::INT_VECTOR};
};

template <> struct infer_type_helper<std::vector<Variant>> {
  static constexpr VariantType value{VariantType::VECTOR};
};

template <> struct infer_type_helper<ObjectId> {
  static constexpr VariantType value{VariantType::OBJECTID};
};
}

template <typename T> constexpr VariantType infer_type() {
  return detail::infer_type_helper<T>::value;
}

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
