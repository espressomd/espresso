/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "Timer.hpp"

#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

namespace Utils {
namespace Timing {

std::vector<std::map<std::string, Utils::Timing::Timer::Stats>>
Timer::get_all_stats() {
  boost::mpi::communicator const &comm = Communication::mpiCallbacks().comm();
  assert(comm.rank() == 0);

  std::vector<std::map<std::string, Utils::Timing::Timer::Stats>> ret(
      comm.size());

  ret[0] = get_stats();

  Communication::mpiCallbacks().call(Timer::mpi_callback);

  for (unsigned i = 1; i < ret.size(); ++i) {
    comm.recv(i, 42, ret[i]);
  }

  return ret;
}

void Timer::mpi_callback(int, int) {
  Communication::mpiCallbacks().comm().send(0, 42, get_stats());
}

std::unordered_map<std::string, Timer> Timer::m_timers;

const Communication::CallbackAdder Timer::cb_adder{Timer::mpi_callback};
}
}
