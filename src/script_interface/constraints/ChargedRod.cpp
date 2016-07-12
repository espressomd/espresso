#include "ChargedRod.hpp"

namespace ScriptInterface {
namespace Constraints {

VariantMap ChargedRod::get_parameters() const {
  return VariantMap{{"center", m_rod.center}, {"lambda", m_rod.lambda}};
}

void ChargedRod::set_parameter(const std::string &name, const Variant &value) {
  std::cout << Communication::mpiCallbacks().comm().rank() << ": "
            << __PRETTY_FUNCTION__ << ": name = " << name << std::endl;

  SET_PARAMETER_HELPER("center", m_rod.center);
  SET_PARAMETER_HELPER("lambda", m_rod.lambda);
}

ParameterMap ChargedRod::all_parameters() const {
  return ParameterMap{
      {"center", Parameter(ParameterType::VECTOR2D, true)},
      {"lambda", Parameter(ParameterType::DOUBLE, true)}};
}
}
}
