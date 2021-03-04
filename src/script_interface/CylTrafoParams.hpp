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

#ifndef SCRIPT_INTERFACE_CYL_TRAFO_PARAMS_HPP
#define SCRIPT_INTERFACE_CYL_TRAFO_PARAMS_HPP

#include "script_interface/ScriptInterface.hpp"

#include "utils/math/cyl_trafo_params.hpp"

namespace ScriptInterface {

class CylTrafoParams : public AutoParameters<CylTrafoParams> {
public:
  CylTrafoParams() {
    add_parameters({{"center", AutoParameter::read_only,
                     [this]() { return m_cyl_trafo_params->center(); }},
                    {"axis", AutoParameter::read_only,
                     [this]() { return m_cyl_trafo_params->axis(); }},
                    {"orientation", AutoParameter::read_only,
                     [this]() { return m_cyl_trafo_params->orientation(); }}});
  }
  std::shared_ptr<::Utils::CylTrafoParams> cyl_trafo_params() {
    return m_cyl_trafo_params;
  }
  void do_construct(VariantMap const &params) override {
    m_cyl_trafo_params = std::make_shared<Utils::CylTrafoParams>(
        get_value_or<Utils::Vector3d>(params, "center",
                                      Utils::Vector3d{{0, 0, 0}}),
        get_value_or<Utils::Vector3d>(params, "axis",
                                      Utils::Vector3d{{0, 0, 1}}),
        get_value_or<Utils::Vector3d>(params, "orientation",
                                      Utils::Vector3d{{1, 0, 0}}));
  }

private:
  std::shared_ptr<Utils::CylTrafoParams> m_cyl_trafo_params;
};
} // namespace ScriptInterface
#endif
