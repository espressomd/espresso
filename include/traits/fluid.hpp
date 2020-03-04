#ifndef OBSERVABLES_FLUID_HPP
#define OBSERVABLES_FLUID_HPP

#include <fluid/fluid.hpp>

namespace Traits {
namespace Fluid {

struct Velocity {
  auto operator()(::Fluid::FluidNode const &fluid) { return fluid.v; }
};

} // namespace Fluid
} // namespace Traits

#endif // OBSERVABLES_FLUID_HPP
