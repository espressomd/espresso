#include "initialize.hpp"

#include "EKSpecies.hpp"

namespace ScriptInterface {
namespace walberla {

void initialize(Utils::Factory<ObjectHandle> *om) {
  om->register_new<EKSpecies>("walberla::EKSpecies");
}

} // namespace walberla
} // namespace ScriptInterface
