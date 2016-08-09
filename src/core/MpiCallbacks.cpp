/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#include <exception>

#include <boost/mpi.hpp>

#include "MpiCallbacks.hpp"

namespace Communication {

void MpiCallbacks::call(int id, int par1, int par2) const {
  std::cout << Communication::mpiCallbacks().comm().rank() << ": "
            << __PRETTY_FUNCTION__ << " id = " << id << std::endl;

  /** Can only be call from master */
  assert(m_comm.rank() == 0);

  /** Check if callback exists */
  if (m_callbacks.find(id) == m_callbacks.end()) {
    throw std::out_of_range("Callback does not exists.");
  }

  int request[3]{id, par1, par2};
  /** Send request to slaves */
  boost::mpi::broadcast(m_comm, request, 3, 0);
}

void MpiCallbacks::call(func_ptr_type fp, int par1, int par2) const {
  /** If the function pointer is invalid, map.at will throw
      an out_of_range exception. */
  const int id = m_func_ptr_to_id.at(fp);

  call(id, par1, par2);
}

int MpiCallbacks::add(const function_type &f) {
  assert(f != nullptr);
  const int id = m_callbacks.add(f);

  assert(m_callbacks.find(id) != m_callbacks.end());

  return id;
}

int MpiCallbacks::add(func_ptr_type fp) {
  assert(fp != nullptr);

  function_type f = function_type(fp);
  const int id = add(f);
  m_func_ptr_to_id[fp] = id;

  return id;
}

void MpiCallbacks::remove(const int id) { m_callbacks.remove(id); }

void MpiCallbacks::slave(int id, int par1, int par2) const {
  std::cout << Communication::mpiCallbacks().comm().rank() << ": "
            << __PRETTY_FUNCTION__ << " id = " << id << std::endl;

  try {
    m_callbacks[id](par1, par2);
  } catch (std::exception &e) {
    std::cout << Communication::mpiCallbacks().comm().rank() << ": "
              << __PRETTY_FUNCTION__ << " id = " << id
              << " failed, what(): " << e.what() << ", aborting... "
              << std::endl;

    mpiCallbacks().comm().abort(253);
  }
}

void MpiCallbacks::abort_loop() const { call(LOOP_ABORT, 0, 0); }

void MpiCallbacks::loop() const {
  for (;;) {
    int request[3];
    /** Communicate callback id and parameters */
    boost::mpi::broadcast(m_comm, request, 3, 0);
    /** id == 0 is loop_abort. */
    if (request[0] == LOOP_ABORT) {
      break;
    } else {
      /** Call the callback */
      slave(request[0], request[1], request[2]);
    }
  }
}

/* We use a singelton callback class for now. */
MpiCallbacks &mpiCallbacks() {
  static boost::mpi::communicator world;
  static MpiCallbacks *m_global_callback = new MpiCallbacks(world);

  return *m_global_callback;
}

} /* namespace Communication */
