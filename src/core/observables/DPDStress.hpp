/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_DPD

#include "Observable.hpp"
#include "dpd.hpp"
#include "system/System.hpp"

#include <cstddef>
#include <vector>

namespace Observables {

class DPDStress : public Observable {
public:
  std::vector<std::size_t> shape() const override { return {3, 3}; }
  std::vector<double>
  operator()(boost::mpi::communicator const &comm) const override {
    return dpd_stress(System::get_system(), comm);
  }
};

} // Namespace Observables

#endif // ESPRESSO_DPD
