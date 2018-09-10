/*
  Copyright (C) 2015-2018 The ESPResSo project

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

#ifndef SCRIPT_INTERFACE_PARAMETERS_HPP
#define SCRIPT_INTERFACE_PARAMETERS_HPP

#include "Variant.hpp"

namespace ScriptInterface {

/**
 * Possible parameter type. Corresponds to the template parameters of @type
 * Variant.
 */
using ParameterType = VariantType;

/**
 * @brief Script interface parameter.
 *
 * Description of a Parameter that can be accessed from
 * the script side.
 */
class Parameter {
public:
  enum { ARBITRARY_LENGTH = 0 };
  Parameter()
      : m_type(ParameterType::BOOL), m_n_elements(0), m_required(false) {}

  Parameter(ParameterType type, bool required)
      : m_type(type), m_required(required) {
    m_n_elements = ARBITRARY_LENGTH;
  }

  Parameter(ParameterType type, size_t n_elements, bool required)
      : m_type(type), m_n_elements(n_elements), m_required(required) {}

  /**
   * The required type of the parameter.
   */
  ParameterType type() const { return m_type; }

  /**
   * The required number of elements of this parameter,
   * if it is a vector type. Otherwise it is ignored.
   */
  size_t n_elements() const { return m_n_elements; }

  /**
   * Whether the parameter is required.
   */
  bool required() const { return m_required; }

private:
  ParameterType m_type;
  size_t m_n_elements;
  bool m_required;
};

} /* namespace ScriptInterface */

#endif
