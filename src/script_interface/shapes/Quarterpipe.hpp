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

#ifndef ESPRESSO_QUARTERPIPE_HPP
#define ESPRESSO_QUARTERPIPE_HPP
#include "Shape.hpp"
#include <shapes/Quarterpipe.hpp>

namespace ScriptInterface {
namespace Shapes {

class Quarterpipe : public Shape {
public:
  Quarterpipe() : m_quarterpipe(new ::Shapes::Quarterpipe()) {
    add_parameters({{"center",
                     [this](Variant const &v) {
                       m_quarterpipe->set_center(get_value<Utils::Vector3d>(v));
                     },
                     [this]() { return m_quarterpipe->center(); }},

                    {"axis",
                     [this](Variant const &v) {
                       m_quarterpipe->set_axis(get_value<Utils::Vector3d>(v));
                     },
                     [this]() { return m_quarterpipe->axis(); }},

                    {"orientation",
                     [this](Variant const &v) {
                       m_quarterpipe->set_orientation(
                           get_value<Utils::Vector3d>(v));
                     },
                     [this]() { return m_quarterpipe->orientation(); }},

                    {"radius",
                     [this](Variant const &v) {
                       m_quarterpipe->set_radius(get_value<double>(v));
                     },
                     [this]() { return m_quarterpipe->radius(); }},

                    {"height",
                     [this](Variant const &v) {
                       m_quarterpipe->set_height(get_value<double>(v));
                     },
                     [this]() { return m_quarterpipe->height(); }}});
  }

  std::shared_ptr<::Shapes::Shape> shape() const override {
    return m_quarterpipe;
  }

private:
  std::shared_ptr<::Shapes::Quarterpipe> m_quarterpipe;
};

} /* namespace Shapes */
} /* namespace ScriptInterface */

#endif
