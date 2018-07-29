#ifndef CORE_CONSTRAINTS_DETAIL_BIND_COUPLING_HPP
#define CORE_CONSTRAINTS_DETAIL_BIND_COUPLING_HPP

namespace FieldCoupling {
namespace detail {
/**
 * @brief Helper class that binds a coupling to a particles
           and propagates Coupling::is_linear.
*/
template <typename Coupling, typename Particle> struct BindCoupling {
  static const constexpr bool is_linear = Coupling::is_linear;

  const Coupling &c;
  const Particle &p;

  BindCoupling(const Coupling &c, const Particle &p) : c(c), p(p) {}

  template <typename T> auto operator()(const T &x) const -> decltype(c(p, x)) {
    return c(p, x);
  }
};

template <typename Coupling, typename Particle>
BindCoupling<Coupling, Particle> make_bind_coupling(const Coupling &c,
                                                    const Particle &p) {
  return BindCoupling<Coupling, Particle>{c, p};
}
}
}

#endif
