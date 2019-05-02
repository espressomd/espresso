/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CORE_CONSTRAINTS_DETAIL_VISCOUS_HPP
#define CORE_CONSTRAINTS_DETAIL_VISCOUS_HPP

#include <utils/Vector.hpp>

namespace FieldCoupling {
namespace Coupling {
class Viscous {
  double m_gamma;

public:
  static constexpr bool is_linear = true;

  Viscous(double gamma) : m_gamma(gamma) {}
  double &gamma() { return m_gamma; }
  double const &gamma() const { return m_gamma; }

  template <typename Particle>
  Utils::Vector3d operator()(Particle const &p,
                             Utils::Vector3d const &field) const {
    return m_gamma * (field - p.m.v);
  }
};
} // namespace Coupling
} // namespace FieldCoupling

#endif
