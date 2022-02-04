#ifndef SCRIPT_INTERFACE_WALBERLA_EK_REACTIONS_HPP
#define SCRIPT_INTERFACE_WALBERLA_EK_REACTIONS_HPP

#include "EKReaction.hpp"
#include "core/grid_based_algorithms/ek_container.hpp"

#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface::walberla {

class EKReactions : public ObjectList<EKReaction> {
  void add_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {}
  void remove_in_core(std::shared_ptr<EKReaction> const &obj_ptr) override {}
};
} // namespace ScriptInterface::walberla

#endif // SCRIPT_INTERFACE_WALBERLA_EK_REACTIONS_HPP
