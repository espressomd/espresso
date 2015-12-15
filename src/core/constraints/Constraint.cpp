#include "Constraint.hpp"
#include "utils/Factory.hpp"

#include "ChargedPlate.hpp"
#include "ChargedRod.hpp"
#include "ExternalMagneticField.hpp"

namespace Constraints {

void initialize_factory() {
  Factory::Instance().register_new("plate", Factory::builder<ChargedPlate>);
  Factory::Instance().register_new("rod", Factory::builder<ChargedRod>);
  Factory::Instance().register_new("ext_magn_field", Factory::builder<ExternalMagneticField>);
}

}
