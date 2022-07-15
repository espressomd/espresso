/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_SHAPES_SHAPE_HPP
#define SCRIPT_INTERFACE_SHAPES_SHAPE_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include <shapes/Shape.hpp>

#include <utils/Vector.hpp>

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Shapes {

class Shape : public AutoParameters<Shape> {
public:
  /**
   * @brief Return the Shape that we are wrapping.
   */
  virtual std::shared_ptr<::Shapes::Shape> shape() const = 0;

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "calc_distance") {
      auto const pos = get_value<Utils::Vector3d>(params.at("position"));
      double dist;
      Utils::Vector3d vec;
      shape()->calculate_dist(pos, dist, vec);
      return std::vector<Variant>{dist, vec};
    }

    if (name == "is_inside") {
      auto const pos = get_value<Utils::Vector3d>(params.at("position"));
      auto is_in = shape()->is_inside(pos);
      return {is_in};
    }

    if (name == "rasterize") {
      auto const grid_size = get_value<Utils::Vector3i>(params.at("grid_size"));
      auto const grid_spacing = get_value<double>(params.at("grid_spacing"));
      auto const grid_offset = get_value<double>(params.at("grid_offset"));
      auto raster = shape()->rasterize(grid_size, grid_spacing, grid_offset);
      return {raster};
    }

    return {};
  }
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
