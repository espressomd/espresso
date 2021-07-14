#include "config.hpp"

#include "EKBoundary.hpp"

#include "core/grid_based_algorithms/ek_boundaries.hpp"
#include "script_interface/ObjectList.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <memory>

#ifndef ESPRESSO_EKBOUNDARIES_HPP
#define ESPRESSO_EKBOUNDARIES_HPP

namespace ScriptInterface {
namespace EKBoundaries {
class EKBoundaries : public ObjectList<EKBoundary> {
  void add_in_core(std::shared_ptr<EKBoundary> const &obj_ptr) override {
    ::EKBoundaries::add(obj_ptr->ekboundary());
  }

  void remove_in_core(std::shared_ptr<EKBoundary> const &obj_ptr) override {
    ::EKBoundaries::remove(obj_ptr->ekboundary());
  }
};
} // namespace EKBoundaries
} /* namespace ScriptInterface */

#endif // ESPRESSO_EKBOUNDARIES_HPP
