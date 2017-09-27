/*
  Copyright (C) 2016,2017 The ESPResSo project

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

#ifndef OBSERVABLES_PIDOBSERVABLE_HPP
#define OBSERVABLES_PIDOBSERVABLE_HPP

#include <vector>
#include "Observable.hpp"

namespace Observables {

// Observable which acts on a given list of particle ids()
class PidObservable : public Observable {
  std::vector<int> m_ids;

public:
  std::vector<int> &ids() { return m_ids; }
  std::vector<int> const &ids() const { return m_ids; }
  int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif
