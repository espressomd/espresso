#ifndef SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP
#define SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP

namespace Interactions {

template <class T> struct LennardJones {
  using value_type = T;
  value_type epsilon;
  value_type sigma;
  value_type energy(value_type r) const {
    return 4.0 * epsilon *
           (std::pow(sigma / r, 12.0) - std::pow(sigma / r, 6.0));
  }
};

} // namespace Interactions

#endif // SRC_INTERACTIONS_INCLUDE_INTERACTIONS_LENNARD_JONES_HPP
