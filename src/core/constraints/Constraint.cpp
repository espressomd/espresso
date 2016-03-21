#include "Constraint.hpp"

#include "ChargedPlate.hpp"
#include "ChargedRod.hpp"
#include "ExternalMagneticField.hpp"
#include "GeometryConstraint.hpp"
#include "ConstraintList.hpp"

namespace Constraints {

using Utils::ParallelFactory;

void initialize_factory() {
  Factory::register_new("plate", Factory::builder<ChargedPlate>);
  Factory::register_new("rod", Factory::builder<ChargedRod>);
  Factory::register_new("ext_magn_field", Factory::builder<ExternalMagneticField>);
  Factory::register_new("geometry_constraint", Factory::builder<GeometryConstraint>);

  ParallelFactory<ConstraintList>::register_new("ConstraintList", ParallelFactory<ConstraintList>::builder<ConstraintList>);
}

}







