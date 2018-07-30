#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP

#include "LBBoundary.hpp"
#include "ScriptInterface.hpp"
#include "ScriptObjectRegistry.hpp"
#include "core/lbboundaries.hpp"

namespace ScriptInterface {
namespace LBBoundaries {
class LBBoundaries : public ScriptObjectRegistry<LBBoundary> {
  void add_in_core(std::shared_ptr<LBBoundary> obj_ptr) override {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    ::LBBoundaries::add(obj_ptr->lbboundary());
#endif
  }

  void remove_in_core(std::shared_ptr<LBBoundary> obj_ptr) override {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    ::LBBoundaries::remove(obj_ptr->lbboundary());
#endif
  }
};
} /* namespace LBBoundaries */
} /* namespace ScriptInterface */

#endif
