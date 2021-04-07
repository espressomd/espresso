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

#ifndef ESPRESSO_HOLLOW_CONICAL_FRUSTUM_HPP
#define ESPRESSO_HOLLOW_CONICAL_FRUSTUM_HPP
#include "Shape.hpp"
#include <script_interface/CylindricalTransformationParameters.hpp>
#include <shapes/HollowConicalFrustum.hpp>

namespace ScriptInterface {
namespace Shapes {

class HollowConicalFrustum : public Shape {
public:
  HollowConicalFrustum() {
    add_parameters(
        {{"cyl_transform_params", m_cyl_transform_params},
         {"r1",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_r1(get_value<double>(v));
          },
          [this]() { return m_hollow_conical_frustum->radius1(); }},
         {"r2",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_r2(get_value<double>(v));
          },
          [this]() { return m_hollow_conical_frustum->radius2(); }},
         {"length",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_length(get_value<double>(v));
          },
          [this]() { return m_hollow_conical_frustum->length(); }},
         {"thickness",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_thickness(get_value<double>(v));
          },
          [this]() { return m_hollow_conical_frustum->thickness(); }},
         {"direction",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_direction(get_value<int>(v));
          },
          [this]() { return m_hollow_conical_frustum->direction(); }},
         {"central_angle",
          [this](Variant const &v) {
            m_hollow_conical_frustum->set_central_angle(get_value<double>(v));
          },
          [this]() { return m_hollow_conical_frustum->central_angle(); }}});
  }

  void do_construct(VariantMap const &params) override {
    set_from_args(m_cyl_transform_params, params, "cyl_transform_params");

    if (m_cyl_transform_params)
      m_hollow_conical_frustum =
          std::make_shared<::Shapes::HollowConicalFrustum>(
              get_value<double>(params, "r1"), get_value<double>(params, "r2"),
              get_value<double>(params, "length"),
              get_value_or<double>(params, "thickness", 0.),
              get_value_or<int>(params, "direction", 1),
              get_value_or<double>(params, "central_angle", 0.),
              m_cyl_transform_params->cyl_transform_params()

          );
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_hollow_conical_frustum;
  }

private:
  std::shared_ptr<::Shapes::HollowConicalFrustum> m_hollow_conical_frustum;
  std::shared_ptr<CylindricalTransformationParameters> m_cyl_transform_params;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
