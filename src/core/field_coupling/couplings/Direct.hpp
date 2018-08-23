#ifndef EXTERNAL_FIELD_COUPLING_DIRECT_HPP
#define EXTERNAL_FIELD_COUPLING_DIRECT_HPP

namespace FieldCoupling {
namespace Coupling {
class Direct {
public:
  static constexpr bool is_linear = true;
  template <typename T, typename Particle>
  const T &operator()(const Particle &, const T &x) const {
    return x;
  }
};
}
}

#endif
