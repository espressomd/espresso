/*
Copyright (C) 2010-2020 The ESPResSo project

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
#pragma once

inline double shear_mode_relaxation_rate(double viscosity) {
  return 2 / (6 * viscosity + 1);
}

inline double odd_mode_relaxation_rate(double shear_relaxation,
                                       double magic_number = 3. / 16.) {
  return (4 - 2 * shear_relaxation) /
         (4 * magic_number * shear_relaxation + 2 - shear_relaxation);
}

inline double viscosity_from_shear_relaxation_rate(double shear_relaxation) {
  return (2 - shear_relaxation) / (6 * shear_relaxation);
}
