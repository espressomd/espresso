#ifndef __MPI_CALLBACKS
#define __MPI_CALLBACKS

#include <functional>
#include <boost/mpi/communicator.hpp>

#include "utils/NumeratedContainer.hpp"

/** Forward declarations so that we don't have to
 * include communication.hpp everywhere,
 * which would slow down compilation considerably.
 */
/*@{*/
void mpi_init(int *argc, char ***argv);
/*@}*/

class MpiCallbacks {
 public:
  /** @brief Function type of callbacks. */
  typedef void (*func_ptr_type)(int, int);
  /** @brief Type of the callback functions. */
  typedef std::function<void (int, int)> function_type;

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
   * @brief Add a new callback.
   *
   * Add a new callback to the system. This is a collective
   * function that must be run on all nodes.
   *
   * @param fp Pointer to the callback function to add.
   * @return An integer id with which the callback can be called.
   **/
  static int add(func_ptr_type fp);
    
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
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  static void call(int id, int par1, int par2);

  /**
   * @brief call a callback.
   *
   * Call the callback id.
   * The method can only be called the master
   * and has the prerequisite that the other nodes are
   * in the MPI loop.
   *
   * @param id The callback to call.
   * @param par1 First parameter to pass to the callback.
   * @param par2 Second parameter to pass to the callback.
   */
  static void call(func_ptr_type fp, int par1, int par2);
    
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
   * @param par2 First parameter to pass to the callback function.
   * @param par2 Second parameter to pass to the callback function.
   */
  static void slave(int id, int par1, int par2);
  
  /**
   * Internal storage for the callback functions.
   */
  static Utils::NumeratedContainer<function_type> m_callbacks;
  /** Mapping of function pointers to ids, so callbabcks can be
   *  called by their point for backward compapbility.
   */
  static std::map<func_ptr_type, int> m_func_ptr_to_id;
};

#endif
