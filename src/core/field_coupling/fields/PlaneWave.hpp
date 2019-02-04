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

#include "gradient_type.hpp"
#include "utils/Vector.hpp"

namespace FieldCoupling {
namespace Fields {
/**
 * @brief A plane wave.
 */
template <typename T, size_t codim> class PlaneWave {
public:
  using value_type = typename decay_to_scalar<Vector<codim, T>>::type;
  using gradient_type = detail::gradient_type<T, codim>;

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

  value_type operator()(const Vector3d &x, T t = 0.) const {
    return m_amplitude * sin(m_k * x - m_omega * t + m_phase);
  }

  gradient_type gradient(const Vector3d &x, T t = 0.) const {
    using Utils::tensor_product;

    return tensor_product(m_amplitude, m_k) *
           cos(m_k * x - m_omega * t + m_phase);
  }

  bool fits_in_box(const Vector3d &) const { return true; }
};
} // namespace Fields
} // namespace FieldCoupling

#endif
