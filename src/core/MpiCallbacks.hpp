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

#ifndef COMMUNICATION_MPI_CALLBACKS
#define COMMUNICATION_MPI_CALLBACKS

#include <functional>

#include <boost/mpi/communicator.hpp>

#include "utils/NumeratedContainer.hpp"

namespace Communication {

/**
 * @brief  The interface of the MPI callback mechanism.
 */
class MpiCallbacks {
public:
  /** Function type of static callbacks. */
  typedef void (*func_ptr_type)(int, int);
  /** Type of the callback functions. */
  typedef std::function<void(int, int)> function_type;

  explicit MpiCallbacks(boost::mpi::communicator const& comm) : m_comm(comm) {
    /** Add a dummy at id 0 for loop abort. */
    m_callbacks.add(function_type());
  }

  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param f The callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  int add(const function_type &f);

  /**
   * @brief Add a new static callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the static callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  int add(func_ptr_type fp);

  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
   *
   * @param id Identifier of the calback to remove.
   */
  void remove(const int id);

  /**
   * @brief call a callback.
   *
   * Call the callback id.
   * The method can only be called on the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  void call(int id, int par1, int par2) const;

  /**
   * @brief call a callback.
   *
   * Call a static callback by pointer.
   * The method can only be called the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id Pointer to static callback (e.g. the function name) to call.
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  void call(func_ptr_type fp, int par1, int par2) const;

  /**
   * @brief Mpi slave loop.
   *
   * This is the callback loop for the slaves. They block
   * on the MPI call and wait for a new callback request
   * coming from the master.
   * This should be run on the slaves and must be running
   * so thet the master can issue call().
   */
  void loop() const;

  /**
   * @brief Abort mpi loop.
   *
   * Make the slaves exit the MPI loop.
   */
  void abort_loop() const;

  /**
   * @brief The boost mpi communicator used by this instance
   */
  boost::mpi::communicator const& comm() const { return m_comm; }

private:
  /**
   * @brief Id for the loop_abort. Has to be 0.
   */
  enum { LOOP_ABORT = 0 };
  /**
   * @brief Callback to integrate with the old callback mechanism.
   *
   * MPI callback for integration with the old callback mechanism,
   * This is called on the slaves by the MPI loop. This should
   * be removed when the new mechanism is used everywhere
   *
   * @param id The id of the callback to run.
   * @param par2 First parameter to pass to the callback function.
   * @param par2 Second parameter to pass to the callback function.
   */
  void slave(int id, int par1, int par2) const;

  /**
   * The MPI communicator used for the callbacks.
   */
  boost::mpi::communicator const& m_comm;

  /**
   * Internal storage for the callback functions.
   */
  Utils::NumeratedContainer<function_type> m_callbacks;

  /**
   * Mapping of function pointers to ids, so static callbacks can be
   * called by their pointer for backward compability.
   */
  std::unordered_map<func_ptr_type, int> m_func_ptr_to_id;
};

/**
 * @brief Initialize the callback singelton.
 * This sets the communicator to use.
 */
  void initialize_callbacks(boost::mpi::communicator const& comm);

/**
 * @brief Returns a reference to the global callback class instance.
 *
 */
MpiCallbacks &mpiCallbacks();

} /* namespace Communication */

#endif
