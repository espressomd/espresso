/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "MpiCallbacks.hpp"

#include <stdexcept>

namespace Communication {
void MpiCallbacks::call(int id, int par1, int par2) const {
  /** Can only be call from master */
  if (m_comm.rank() != 0) {
    throw std::logic_error("Callbacks can only be invoked on rank 0.");
  }

  /** Check if callback exists */
  if (m_callbacks.find(id) == m_callbacks.end()) {
    throw std::out_of_range("Callback does not exists.");
  }

  /** Send request to slaves */
  boost::mpi::broadcast(m_comm, id, 0);
  detail::tuple<int,int> params{{par1, par2}};
  boost::mpi::broadcast(m_comm, params, 0);
}

void MpiCallbacks::call(func_ptr_type fp, int par1, int par2) const {
  /** If the function pointer is invalid, map.at will throw
      an out_of_range exception. */
  const int id = m_func_ptr_to_id.at(reinterpret_cast<void *>(fp));

  call(id, par1, par2);
}

int MpiCallbacks::add(const function_type &f) {
  const int id = m_callbacks.add(detail::make_model(f));

  assert(m_callbacks.find(id) != m_callbacks.end());

  return id;
}

int MpiCallbacks::add(func_ptr_type fp) {
  assert(fp);

  const int id = m_callbacks.add(detail::make_model(fp));
  m_func_ptr_to_id[reinterpret_cast<void*>(fp)] = id;

  return id;
}

void MpiCallbacks::remove(const int id) { m_callbacks.remove(id); }

void MpiCallbacks::abort_loop() const { call(LOOP_ABORT, 0, 0); }

void MpiCallbacks::loop() const {
  for (;;) {
    int request;
    /** Communicate callback id and parameters */
    boost::mpi::broadcast(m_comm, request, 0);
    /** id == 0 is loop_abort. */
    if (request == LOOP_ABORT) {
      break;
    } else {
      /** Call the callback */
      m_callbacks[request]->operator()(m_comm);
    }
  }
}
} /* namespace Communication */
