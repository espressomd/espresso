#ifndef SCRIPT_INTERFACE_WALBERLA_INITIALIZE_HPP
#define SCRIPT_INTERFACE_WALBERLA_INITIALIZE_HPP

#include "script_interface/ObjectHandle.hpp"
#include "utils/Factory.hpp"

namespace ScriptInterface {
namespace walberla {
void initialize(Utils::Factory<ObjectHandle> *om);
} // namespace walberla
} // namespace ScriptInterface

#endif // SCRIPT_INTERFACE_WALBERLA_INITIALIZE_HPP
