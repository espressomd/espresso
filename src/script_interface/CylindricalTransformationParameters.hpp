/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_CYL_TRANSFORM_PARAMS_HPP
#define SCRIPT_INTERFACE_CYL_TRANSFORM_PARAMS_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include <utils/math/cylindrical_transformation_parameters.hpp>

#include <memory>
#include <stdexcept>

namespace ScriptInterface {

class CylindricalTransformationParameters
    : public AutoParameters<CylindricalTransformationParameters> {
public:
  CylindricalTransformationParameters() {
    add_parameters({{"center", AutoParameter::read_only,
                     [this]() { return m_transform_params->center(); }},
                    {"axis", AutoParameter::read_only,
                     [this]() { return m_transform_params->axis(); }},
                    {"orientation", AutoParameter::read_only,
                     [this]() { return m_transform_params->orientation(); }}});
  }
  std::shared_ptr<::Utils::CylindricalTransformationParameters>
  cyl_transform_params() const {
    return m_transform_params;
  }
  void do_construct(VariantMap const &params) override {
    auto n_params = params.size();
    switch (n_params) {
    case 0:
      m_transform_params =
          std::make_shared<Utils::CylindricalTransformationParameters>();
      break;
    case 2:
      m_transform_params =
          std::make_shared<Utils::CylindricalTransformationParameters>(
              get_value<Utils::Vector3d>(params, "center"),
              get_value<Utils::Vector3d>(params, "axis"));
      break;
    case 3:
      m_transform_params =
          std::make_shared<Utils::CylindricalTransformationParameters>(
              get_value<Utils::Vector3d>(params, "center"),
              get_value<Utils::Vector3d>(params, "axis"),
              get_value<Utils::Vector3d>(params, "orientation"));
      break;
    default:
      throw std::runtime_error("Provide either no arguments, center and axis, "
                               "or center and axis and orientation");
    }
  }

private:
  std::shared_ptr<Utils::CylindricalTransformationParameters>
      m_transform_params;
};
} // namespace ScriptInterface
#endif
