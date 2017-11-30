#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP

#include "LBBoundary.hpp"
#include "ScriptInterface.hpp"
#include "ScriptObjectRegistry.hpp"
#include "config.hpp"
#include "core/lbboundaries.hpp"

namespace ScriptInterface {
namespace LBBoundaries {

class LBBoundaries : public ScriptObjectRegistry<LBBoundary> {
  virtual void add_in_core(std::shared_ptr<LBBoundary> obj_ptr) override {
    ::LBBoundaries::lbboundaries.push_back(obj_ptr->lbboundary());
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
    ::LBBoundaries::lb_init_boundaries();
#endif
  }

  virtual void remove_in_core(std::shared_ptr<LBBoundary> obj_ptr) override {
    auto it = std::find(std::begin(::LBBoundaries::lbboundaries),
                        std::end(::LBBoundaries::lbboundaries),
                        obj_ptr->lbboundary());

    if (it != std::end(::LBBoundaries::lbboundaries)) {
      ::LBBoundaries::lbboundaries.erase(it);
    }
  }
};
} /* namespace LBBoundaries */
} /* namespace ScriptInterface */

#endif
