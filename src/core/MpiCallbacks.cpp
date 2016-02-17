#include "MpiCallbacks.hpp"
#include "communication.hpp"

#include <boost/mpi.hpp>

void MpiCallbacks::call(int id, int par) {
  assert((id >= 0) && (id < m_callbacks.size()));
  mpi_call(MpiCallbacks::slave, id, par);
}


int MpiCallbacks::add(const Function &f) {
  assert(f != nullptr);
  const int id = m_callbacks.size();
  m_callbacks.push_back(f);

  return id;
}


void MpiCallbacks::slave(int id, int par) {
  assert(id < m_callbacks.size());
  m_callbacks[id](par);
}

boost::mpi::communicator MpiCallbacks::mpi_comm;
