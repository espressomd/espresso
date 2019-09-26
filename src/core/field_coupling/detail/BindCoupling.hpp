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
#ifndef CORE_CONSTRAINTS_DETAIL_BIND_COUPLING_HPP
#define CORE_CONSTRAINTS_DETAIL_BIND_COUPLING_HPP

namespace FieldCoupling {
namespace detail {
/**
 * @brief Helper class that binds a coupling to a particles
           and propagates Coupling::is_linear.
*/
template <typename Coupling, typename Particle> struct BindCoupling {
  static const constexpr bool is_linear = Coupling::is_linear;

  const Coupling &c;
  const Particle &p;

  BindCoupling(const Coupling &c, const Particle &p) : c(c), p(p) {}

  template <typename T> auto operator()(const T &x) const -> decltype(c(p, x)) {
    return c(p, x);
  }
};

template <typename Coupling, typename Particle>
BindCoupling<Coupling, Particle> make_bind_coupling(const Coupling &c,
                                                    const Particle &p) {
  return BindCoupling<Coupling, Particle>{c, p};
}
} // namespace detail
} // namespace FieldCoupling

#endif
