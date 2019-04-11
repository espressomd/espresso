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
#ifndef CORE_CONSTRAINTS_DETAIL_SCALED_HPP
#define CORE_CONSTRAINTS_DETAIL_SCALED_HPP

#include <unordered_map>

namespace FieldCoupling {
namespace Coupling {
class Scaled {
  std::unordered_map<int, double> m_scales;
  double m_default;

public:
  static constexpr bool is_linear = true;

  template <typename ScalesRef>
  Scaled(ScalesRef &&scales, double default_val)
      : m_scales(std::forward<ScalesRef>(scales)), m_default(default_val) {}

  double &default_scale() { return m_default; }
  double const &default_scale() const { return m_default; }
  std::unordered_map<int, double> &particle_scales() { return m_scales; }
  std::unordered_map<int, double> const &particle_scales() const {
    return m_scales;
  }

private:
  template <typename Particle> double scale(Particle const &p) const {
    auto const &val = m_scales.find(p.identity());
    if (val != m_scales.end())
      return val->second;

    return m_default;
  }

public:
  template <typename T, typename Particle>
  T operator()(const Particle &p, const T &x) const {
    return scale(p) * x;
  }
};
} // namespace Coupling
} // namespace FieldCoupling

#endif
