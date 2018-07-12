#ifndef CORE_CONSTRAINTS_DETAIL_CHARGE_HPP
#define CORE_CONSTRAINTS_DETAIL_CHARGE_HPP

namespace FieldCoupling {
namespace Coupling {
class Charge {
public:
  static constexpr const bool is_linear = true;

  template <typename T, typename Particle>
  T operator()(const Particle &p, const T &x) const {
    return p.p.q * x;
  }
};
}
}

#endif
