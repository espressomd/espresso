#ifndef SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP
#define SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

namespace Interactions {

template <class T> struct LennardJones {
  using value_type = T;
  value_type epsilon;
  value_type sigma;

  value_type energy(value_type r) const {
    return 4.0 * epsilon *
           (Utils::int_pow<12>(sigma / r) - Utils::int_pow<6>(sigma / r));
  }

  Utils::Vector<value_type, 3>
  force(value_type r, Utils::Vector<value_type, 3> const &r12) const {
    return r12 *
           (24. * epsilon * Utils::int_pow<6>(sigma) *
            (Utils::int_pow<6>(r) - 2. * Utils::int_pow<6>(sigma))) /
           Utils::int_pow<13>(r);
  }
};

} // namespace Interactions

#endif // SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP
