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
#ifndef EXTERNAL_FIELD_DETAIL_BASE_HPP
#define EXTERNAL_FIELD_DETAIL_BASE_HPP

#include <utility>

namespace FieldCoupling {
namespace detail {
template <typename Coupling, typename Field> class Base {
protected:
  Coupling m_coupling;
  Field m_field;

public:
  template <typename CouplingRef, typename FieldRef>
  Base(CouplingRef &&coupling, FieldRef &&field)
      : m_coupling(std::forward<CouplingRef>(coupling)),
        m_field(std::forward<FieldRef>(field)) {}

  Coupling const &coupling() const { return m_coupling; }
  Field const &field() const { return m_field; }
};
} // namespace detail
} // namespace FieldCoupling
#endif
