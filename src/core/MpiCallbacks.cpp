#include "MpiCallbacks.hpp"
#include "communication.hpp"

#include <boost/mpi.hpp>

void MpiCallbacks::call(int id, int par) {  
  /** Can only be call from master, otherwise it would
      be called from a callback, which would lead to MPI deadlock. */
  assert(mpi_comm.rank() == 0);
  mpi_call(MpiCallbacks::slave, id, par);
}

int MpiCallbacks::add(function_type &f) {
  assert(f != nullptr);
  const int id = m_callbacks.add(f);

  return id;
}

void MpiCallbacks::remove(const int id) {
  m_callbacks.remove(id);
}

void MpiCallbacks::slave(int id, int par) {
  m_callbacks[id](par);
}

/** Initialize static members */
boost::mpi::communicator MpiCallbacks::mpi_comm;
ObjectContainer<MpiCallbacks::function_type> MpiCallbacks::m_callbacks;
