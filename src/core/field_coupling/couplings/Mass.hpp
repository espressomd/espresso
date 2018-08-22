#ifndef CORE_CONSTRAINTS_DETAIL_MASS_HPP
#define CORE_CONSTRAINTS_DETAIL_MASS_HPP

namespace FieldCoupling {
namespace Coupling {
class Mass {
public:
  static constexpr const bool is_linear = true;

  template <typename T, typename Particle>
  T operator()(const Particle &p, const T &x) const {
    return p.p.mass * x;
  }
};
}
}

#endif
