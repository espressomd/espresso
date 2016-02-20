#ifndef __MPI_CALLBACKS
#define __MPI_CALLBACKS

#include <functional>
#include <boost/mpi/communicator.hpp>

#include "ObjectContainer.hpp"

/** Forward declarations so that we don't have to
 * include communication.hpp everywhere,
 * which would slow down compilation considerably.
 */
/*@{*/
void mpi_init(int *argc, char ***argv);
/*@}*/

class MpiCallbacks {
 public:
  /** @brief Type of the callback functions */
  typedef std::function<void (int)> function_type;

  /**
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param f The callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  static int add(function_type &f);
  
  /**
   * @brief Remove callback.
   *
   * Remove the callback id from the callback list.
   * This is a collective call that must be run on all node.
   *
   * @param id Identifier of the calback to remove.
   */
  static void remove(const int id);
  
  /**
   * @brief call a callback.
   *
   * Call the callback id.
   * The method can only be called the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param par The parameter to pass to the callback.
   */
  static void call(int id, int par);
  
  /**
   * The MPI communicator used for the callbacks.   
   */
  static boost::mpi::communicator mpi_comm;
  
 private:
  
  /**
   * mpi_init is friended so that it can insert the constant
   * callbacks form communication.cpp at initialization.
   * This should be removed when the dynamic mechanism to add
   * callbacks is used everywhere.
   */
  friend void mpi_init(int *argc, char ***argv);
  
  /**
   * @briefcCallback to integrate with the old callback mechanism.
   *
   * MPI callback for integration with the old callback mechanism,
   * This is called on the slaves by the MPI loop. This should
   * be removed when the new mechanism is used everywhere
   *
   * @param id The id of the callback to run.
   * @param par THe parameter to pass to the callback function.
   */
  static void slave(int id, int par);
  
  /**
   * Internal storage for the callback functions.
   */
  static ObjectContainer<function_type> m_callbacks;
};

#endif
