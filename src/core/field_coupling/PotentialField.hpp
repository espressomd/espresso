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
#ifndef FIELD_COUPLING_EXTERNAL_POTENTIAL_HPP
#define FIELD_COUPLING_EXTERNAL_POTENTIAL_HPP

#include "detail/Base.hpp"
#include "detail/BindCoupling.hpp"

#include <utils/Vector.hpp>

namespace FieldCoupling {
template <typename Coupling, typename Field>
class PotentialField : public detail::Base<Coupling, Field> {
  using Base = detail::Base<Coupling, Field>;

  using Base::m_coupling;
  using Base::m_field;

public:
  using Base::Base;

  template <typename Particle>
  Utils::Vector3d force(const Particle &p, const Utils::Vector3d &folded_pos,
                        double t) const {
    using detail::make_bind_coupling;
    return m_coupling(p, -m_field.jacobian(folded_pos, t));
  }

  template <typename Particle>
  double energy(const Particle &p, const Utils::Vector3d &folded_pos,
                double t) const {
    return m_coupling(p, m_field(folded_pos, t));
  }
};
} /* namespace FieldCoupling */

#endif
