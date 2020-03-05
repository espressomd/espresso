#ifndef INCLUDE_FLUID_FLUID_HPP
#define INCLUDE_FLUID_FLUID_HPP

namespace Fluid {

struct FluidNode {
  double v;
};

} // namespace Fluid

namespace Traits {
namespace Fluid {

struct Velocity {
  auto operator()(::Fluid::FluidNode const &fluid) { return fluid.v; }
};

} // namespace Fluid
} // namespace Traits

#endif // INCLUDE_FLUID_FLUID_HPP
