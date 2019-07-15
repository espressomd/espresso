#ifndef ESPRESSO_SCRIPT_INTERFACE_OBJECTSTATE_HPP
#define ESPRESSO_SCRIPT_INTERFACE_OBJECTSTATE_HPP

#include "Variant.hpp"

#include <string>
#include <utility>
#include <vector>

namespace ScriptInterface {
/**
 * @brief State of an object ready for serialization.
 *
 * This specifies the internal serialization format and
 * should not be used outside of the class.
 */
struct ObjectState {
  std::string name;
  CreationPolicy policy;
  PackedMap params;
  std::vector<std::pair<ObjectId, std::string>> objects;
  std::string internal_state;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &name &policy &params &objects &internal_state;
  }
};
} // namespace ScriptInterface

#endif
