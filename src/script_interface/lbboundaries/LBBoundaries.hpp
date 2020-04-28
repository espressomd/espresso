/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP
#define SCRIPT_INTERFACE_LBBOUNDARIES_LBBOUNDARIES_HPP

#include "LBBoundary.hpp"
#include "core/grid_based_algorithms/lb_boundaries.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/ScriptObjectRegistry.hpp"
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
