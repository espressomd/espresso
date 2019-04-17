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
#ifndef CORE_EXTERNAL_FIELD_FIELDS_PLAIN_WAVE_HPP
#define CORE_EXTERNAL_FIELD_FIELDS_PLAIN_WAVE_HPP

#include "jacobian_type.hpp"
#include "utils/Vector.hpp"

namespace FieldCoupling {
namespace Fields {
/**
 * @brief A plane wave.
 *
 * A time dependent plane wave, with a certain (vector-valued)
 * amplitude, wave vector frequency and phase.
 *
 * See https://en.wikipedia.org/wiki/Plane_wave
 */
template <typename T, size_t codim> class PlaneWave {
public:
  using value_type =
      typename Utils::decay_to_scalar<Utils::Vector<T, codim>>::type;
  using jacobian_type = detail::jacobian_type<T, codim>;

private:
  value_type m_amplitude;
  value_type m_k;
  T m_omega;
  T m_phase;

public:
  PlaneWave(const value_type amplitude, const value_type &k, T omega, T phase)
      : m_amplitude(amplitude), m_k(k), m_omega(omega), m_phase(phase) {}

  value_type &amplitude() { return m_amplitude; }
  value_type &k() { return m_k; }
  T &omega() { return m_omega; }
  T &phase() { return m_phase; }

  /**
   * brief Evaluate field.
   *
   * @param x Where?
   * @param t When?
   * @return Value of the field at point x and time t.
   */
  value_type operator()(const Utils::Vector3d &x, T t = 0.) const {
    return m_amplitude * sin(m_k * x - m_omega * t + m_phase);
  }

  /**
   * brief Evaluate the Jacobian matrix of the field.
   *
   * See https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
   * In the special case of a scalar field, the Jacobian is the gradient of
   * the field.
   *
   * @param x Where?
   * @param t When?
   * @return Jacobian matrix
   */
  jacobian_type jacobian(const Utils::Vector3d &x, T t = 0.) const {
    using Utils::tensor_product;

    return tensor_product(m_amplitude, m_k) *
           cos(m_k * x - m_omega * t + m_phase);
  }

  bool fits_in_box(const Utils::Vector3d &) const { return true; }
};
} // namespace Fields
} // namespace FieldCoupling

#endif
