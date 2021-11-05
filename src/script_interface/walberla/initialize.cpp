#include <script_interface/ObjectHandle.hpp>
#include <utils/Factory.hpp>

#include "FluidWalberla.hpp"
#include "LatticeWalberla.hpp"

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<LatticeWalberla>("walberla::LatticeWalberla");
  om->register_new<FluidWalberla>("walberla::FluidWalberla");
}

} // namespace ScriptInterface::walberla
