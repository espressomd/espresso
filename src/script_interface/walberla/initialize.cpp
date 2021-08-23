#include "initialize.hpp"

#include "EKContainer.hpp"
#include "EKFFT.hpp"
#include "EKNone.hpp"
#include "EKSpecies.hpp"

namespace ScriptInterface::walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<EKContainer>("walberla::EKContainer");
  om->register_new<EKSpecies>("walberla::EKSpecies");
  om->register_new<EKFFT>("walberla::EKFFT");
  om->register_new<EKNone>("walberla::None");
}

} // namespace ScriptInterface::walberla
