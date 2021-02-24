#ifndef SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LINEAR_HPP
#define SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LINEAR_HPP

#include <utils/Vector.hpp>

namespace Interactions {

template <class T> struct Linear {
  using value_type = T;

  value_type m_c1;
  value_type m_c2;

  value_type energy(value_type r) const { return -m_c1 * r + m_c2; }

  Utils::Vector<value_type, 3>
  force(value_type r, Utils::Vector<value_type, 3> const &r12) const {
    return r12 * m_c1;
  }
};

} // namespace Interactions

#endif // SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LINEAR_HPP
