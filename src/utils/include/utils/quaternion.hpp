#ifndef SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
#define SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP

#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/quat_traits.hpp>

#include "utils/Vector.hpp"

namespace Utils {
template <typename T> class Quaternion : public Vector<T, 4> {
  using base = Vector<T, 4>;

public:
  using base::base;
  Quaternion &normalize() {
    boost::qvm::normalize(*this);
    return *this;
  }

  T norm() { return boost::qvm::mag(*this); }
};
} // namespace Utils

namespace boost {

namespace qvm {

template <class T> struct quat_traits<Utils::Quaternion<T>> {
  using quat_type = Utils::Quaternion<T>;
  using scalar_type = typename quat_type::value_type;

  template <std::size_t I>
  static constexpr inline scalar_type &write_element(quat_type &q) {
    return q[I];
  }

  template <std::size_t I>
  static constexpr inline scalar_type read_element(quat_type const &q) {
    return q[I];
  }
};

} // namespace qvm
} // namespace boost
#endif // SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
