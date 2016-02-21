#include "MpiCallbacks.hpp"
#include "communication.hpp"

#include <boost/mpi.hpp>

void MpiCallbacks::call(int id, int par1, int par2) {  
  /** Can only be call from master */
  assert(mpi_comm.rank() == 0);
  int request[3]{id, par1, par2};
  /** Send request to slaves */
  boost::mpi::broadcast(mpi_comm, request, 3, 0);
}

void MpiCallbacks::call(func_ptr_type fp, int par1, int par2) {  
  cont int id = m_func_ptr_to_id[fp];
  call(id, par1, par2);
}

int MpiCallbacks::add(function_type &f) {
  assert(f != nullptr);
  const int id = m_callbacks.add(f);

  return id;
}

int MpiCallbacks::add(func_ptr_type fp) {
  assert(f != nullptr);
  
  const int id = add(std::function<function_type>(fp);

  return id;
}

void MpiCallbacks::remove(const int id) {
  m_callbacks.remove(id);
}

void MpiCallbacks::slave(int id, int par1, int par2) {
  m_callbacks[id](par1, par2);
}

void MpiCallbacks::loop() {
  for(;;) {
    int request[3];
    /** Communicate callback id and parameters */
    boost::mpi::broadcast(mpi_comm, request, 3, 0);
    /** Call the callback */
    call(request[0], request[1], request[2]);
  }
}

/** Initialize static members */
boost::mpi::communicator MpiCallbacks::mpi_comm;
Utils::NumeratedContainer<MpiCallbacks::function_type> MpiCallbacks::m_callbacks;
std::map<func_ptr_type, int> MpiCallbacks::m_func_ptr_to_id;
