#include "initialize.hpp"

#include "EKContainer.hpp"
#include "EKFFT.hpp"
#include "EKSpecies.hpp"

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<EKContainer>("walberla::EKContainer");
  om->register_new<EKSpecies>("walberla::EKSpecies");
  om->register_new<EKFFT>("walberla::EKFFT");
}

} // namespace ScriptInterface::walberla
