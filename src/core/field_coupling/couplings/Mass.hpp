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
#ifndef CORE_CONSTRAINTS_DETAIL_MASS_HPP
#define CORE_CONSTRAINTS_DETAIL_MASS_HPP

namespace FieldCoupling {
namespace Coupling {
class Mass {
public:
  static constexpr const bool is_linear = true;

  template <typename T, typename Particle>
  T operator()(const Particle &p, const T &x) const {
    return p.p.mass * x;
  }
};
} // namespace Coupling
} // namespace FieldCoupling

#endif
