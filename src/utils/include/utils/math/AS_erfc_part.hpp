/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017,2018 The ESPResSo
  project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef UTILS_MATH_AC_ERF_PART_HPP
#define UTILS_MATH_AC_ERF_PART_HPP

namespace Utils {
/** approximates \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
    Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
    (9. ed.), chapter 7 */
template <typename T> constexpr T AS_erfc_part(T d) {
  T const constexpr a1 = 0.254829592;
  T const constexpr a2 = -0.284496736;
  T const constexpr a3 = 1.421413741;
  T const constexpr a4 = -1.453152027;
  T const constexpr a5 = 1.061405429;
  T const constexpr p = 0.3275911;

  auto const t = 1.0 / (1.0 + p * d);

  return t * (a1 + t * (a2 + t * (a3 + t * (a4 + t * a5))));
}
} // namespace Utils

#endif
