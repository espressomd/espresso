#ifndef SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
#define SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP

#include <boost/qvm/deduce_scalar.hpp>
#include <boost/qvm/quat_operations.hpp>
#include <boost/qvm/quat_traits.hpp>

#include "utils/Vector.hpp"

namespace Utils {
template <typename T> class Quaternion : private Vector<T, 4> {
  using base = Vector<T, 4>;
  friend class boost::serialization::access;

public:
  using base::base;
  using base::operator[];
  using base::begin;
  using base::data;
  using base::end;
  using typename base::value_type;
  Quaternion<T>() = default;
  Quaternion<T>(Quaternion<T> const &) = default;
  Quaternion<T> &operator=(Quaternion<T> const &) = default;
  Quaternion<T>(std::initializer_list<T> values) {
    std::copy(values.begin(), values.end(), begin());
  };

  Quaternion &normalize() {
    boost::qvm::normalize(*this);
    return *this;
  }

  T norm() const { return boost::qvm::mag(*this); }
  T norm2() const { return boost::qvm::mag_sqr(*this); }

  static Quaternion<T> identity() { return boost::qvm::identity_quat<T>(); }

  static Quaternion<T> zero() { return boost::qvm::zero_quat<T>(); }
};

template <typename T> T dot(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::dot(a, b);
}

template <typename T>
Quaternion<T> operator*(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator*(a, b);
}

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value, bool> = true>
Quaternion<T> operator*(const Quaternion<T> &a, const U &b) {
  return boost::qvm::operator*(a, b);
}

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value, bool> = true>
Quaternion<T> operator*(const U &b, const Quaternion<T> &a) {
  return boost::qvm::operator*(a, b);
}

template <typename T>
Quaternion<T> operator/(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator/(a, b);
}

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value, bool> = true>
Quaternion<T> operator/(const Quaternion<T> &a, const U &b) {
  return boost::qvm::operator/(a, b);
}

template <typename T, typename U,
          std::enable_if_t<std::is_arithmetic<U>::value, bool> = true>
Quaternion<T> operator/=(Quaternion<T> &a, const U &b) {
  return boost::qvm::operator/=(a, b);
}

// template <typename T>
// Quaternion<T> operator/=(Quaternion<T> &a, const Quaternion<T> &b) {
//  return boost::qvm::operator/=(a, b);
//}

template <typename T>
Quaternion<T> operator+(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator+(a, b);
}

template <typename T>
Quaternion<T> operator+=(Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator+=(a, b);
}

template <typename T>
Quaternion<T> operator-(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator-(a, b);
}

template <typename T>
Quaternion<T> operator-=(Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator-=(a, b);
}

template <typename T>
bool operator==(const Quaternion<T> &a, const Quaternion<T> &b) {
  return boost::qvm::operator==(a, b);
}
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

template <typename T>
struct deduce_scalar<Utils::Quaternion<T>, Utils::Quaternion<T>> {
  using type = T;
};

} // namespace qvm
} // namespace boost
#endif // SRC_UTILS_INCLUDE_UTILS_QUATERNION_HPP
