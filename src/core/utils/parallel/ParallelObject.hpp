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

#ifndef UTILS_PARALLEL_OBJECT_HPP
#define UTILS_PARALLEL_OBJECT_HPP

#include "MpiCallbacks.hpp"

namespace Utils {
template <typename T> /* The type we are wrapping */
class ParallelObject {
  ParallelObject() { Communication::mpiCallbacks().call(m_callback_id, 0, 0); }
  /**
   * Decay into wrapped type.
   */
  operator &T() { return m_d; }

private:
  enum CallbackAction { CREATE };

  friend Communication::MpiCallbacks;
  static void mpi_callback(int a, int b) {
    switch (action) {
    case CREATE:
      /* Create an instance */
      new T();
      break;
    }
  }

  static int m_callback_id;

  T m_d;
};

  template<typename T>
      ParallelObject::m_callback_id = Communication::mpiCallbacks().add(ParallelObject<T>::mpi_callback;

  } /* namespace Utils */

#endif
