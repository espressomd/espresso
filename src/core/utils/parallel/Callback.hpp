/*
  Copyright (C) 2016 The ESPResSo project

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

#ifndef UTILS_PARALLEL_INSTANCE_CALLBACK_HPP
#define UTILS_PARALLEL_INSTANCE_CALLBACK_HPP

#include "core/MpiCallbacks.hpp"

namespace Utils {
namespace Parallel {
/**
 * @brief Add a mpi callback to a class.
 *
 * This is a RAII class to register a mpi callbcack
 * per instance. The callback has the same livetime as
 * the instance: It is created in the constructor
 * and removed in the destructor.
 */
class Callback {
public:
  Callback(Communication::MpiCallbacks &cb,
           const Communication::MpiCallbacks::function_type &callback) : m_cb(cb) {
    m_callback_id = m_cb.add(callback);
  }

  ~Callback() { m_cb.remove(m_callback_id); }

  /**
   * @brief Run the callback function on the slave.
   *
   * The callback is not run on the calling node.
   */
  void call(int a = 0, int b = 0) {
    m_cb.call(m_callback_id, a, b);
  }

private:
  /* Callback system we're on */
  Communication::MpiCallbacks & m_cb;
  /* Id of the encapsulated callback */
  int m_callback_id;
};

} /* namespace Parallel */
} /* namespace Utils */

#endif
