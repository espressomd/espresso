/*
  Copyright (C) 2015,2016 The ESPResSo project

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

#ifndef UTILS_PARALLEL_PARALLEL_OBJECT_HPP
#define UTILS_PARALLEL_PARALLEL_OBJECT_HPP

#include "MpiCallbacks.hpp"

namespace Utils {
namespace Parallel {

template <typename T> /* The type we are wrapping */
class ParallelObject {
public:
  ParallelObject() = delete;

  static void register_callback(Communication::MpiCallbacks & cb) {
    cb.add(&mpi_callback);
  }

  static void make(Communication::MpiCallbacks & cb) {
    cb.call(&mpi_callback);
  }

private:
  /* Supported callback types. Currently we can only create new instances. */
  enum CallbackAction { CREATE };

  friend Communication::MpiCallbacks;
  static void mpi_callback(int action, int) {
    switch (action) {
    case CREATE:
      /* Create an instance */
      new T();
      break;
    }
  }
};

} /* namespace Parallel */
} /* namespace Utils */
#endif
