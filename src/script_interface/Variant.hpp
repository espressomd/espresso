#ifndef SCRIPT_INTERFACE_VARIANT_HPP
#define SCRIPT_INTERFACE_VARIANT_HPP

#include <boost/variant.hpp>

#include "None.hpp"
#include "core/Vector.hpp"
#include "utils/AutoObjectId.hpp"

#include <boost/serialization/map.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>

#include <map>

namespace ScriptInterface {
class ScriptInterfaceBase;
using ObjectId = Utils::ObjectId<ScriptInterfaceBase>;

/**
 * @brief None-"literal".
 */
constexpr const None none{};

/**
 * @brief Possible types for parameters.
 */
typedef boost::make_recursive_variant<
    None, bool, int, double, std::string, std::vector<int>, std::vector<double>,
    ObjectId, std::vector<boost::recursive_variant_>>::type Variant;

/**
 * @brief Human readable names for the types in the variant.
 *
 * This should be in the same order as the corresponding types in the definition
 * of Variant. (e.g. Variant(bool{}).which() == VariantType::BOOL).
 */
enum class VariantType {
  NONE = 0,
  BOOL,
  INT,
  DOUBLE,
  STRING,
  INT_VECTOR,
  DOUBLE_VECTOR,
  OBJECTID,
  VECTOR
};

typedef std::map<std::string, Variant> VariantMap;

namespace detail {
/**
 * @brief Implementation of @f infer_type.
 *
 * Helper struct is needed because particle specialization
 * of functions is not allowed. Every specialization deals
 * with a specific type, to extend just add a new one.
 */
template <typename T> struct infer_type_helper {
  static_assert(sizeof(T) == 0, "Type T is not contained in Variant.");
};

template <> struct infer_type_helper<None> {
  static constexpr VariantType value{VariantType::NONE};
};

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

/* Deduce std::shared_ptr<ScriptInterfaceBase> as ObjectId */
template <typename T> struct infer_type_helper<std::shared_ptr<T>> {
  static_assert(std::is_base_of<ScriptInterfaceBase, T>::value, "");
  static constexpr VariantType value{VariantType::OBJECTID};
};
}

/**
 * @brief Infer the variant type id from the c++ type.
 *
 * infer_type<int>() returns VariantType::INT an so on.
 */
template <typename T> constexpr VariantType infer_type() {
  return detail::infer_type_helper<T>::value;
}

/**
 * @brief Get a string representation of the current type of a variant.
 */
std::string get_type_label(Variant const &);

/**
 * @brief Get a string representation of VariantType.
 */
std::string get_type_label(VariantType);

bool is_none(Variant const &v);
bool is_bool(Variant const &v);
bool is_int(Variant const &v);
bool is_string(Variant const &v);
bool is_double(Variant const &v);
bool is_int_vector(Variant const &v);
bool is_double_vector(Variant const &v);
bool is_objectid(Variant const &v);
bool is_vector(Variant const &v);

/**
 * @brief Combine doubles/ints into respective vectors.
 *
 * This function checks recursively for vectors of variants
 * that are all doubles or ints, and if so combines them into
 * a single double/int vectors. So if you have a variant with
 * types { { INT, INT, INT }, STRING, {DOUBLE, DOUBLE}},
 * this would be reduced to { INT_VECTOR, STRING, DOUBLE_VECTOR }.
 *
 * @param v The variant to transform.
 */
void transform_vectors(Variant &v);

/**
 * @brief Recursivly print the type of a variant.
 */
std::string print_variant_types(Variant const &v);

} /* namespace ScriptInterface */

#endif
