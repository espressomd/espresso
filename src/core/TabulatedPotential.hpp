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
#ifndef CORE_TABULATED_POTENTIAL_HPP
#define CORE_TABULATED_POTENTIAL_HPP

#include "utils/linear_interpolation.hpp"
#include "utils/serialization/List.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <vector>

struct TabulatedPotential {
  double minval = -1.0;
  double maxval = -1.0;
  double invstepsize = 0.0;
  std::vector<double> force_tab;
  std::vector<double> energy_tab;

  double force(double x) const {
    assert(x <= maxval);
    return Utils::linear_interpolation(force_tab, invstepsize, minval, x);
  }

  double energy(double x) const {
    assert(x <= maxval);
    return Utils::linear_interpolation(energy_tab, invstepsize, minval, x);
  }

  double cutoff() const { return maxval; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &minval;
    ar &maxval;
    ar &invstepsize;
    ar &force_tab;
    ar &energy_tab;
  }
};

#endif
