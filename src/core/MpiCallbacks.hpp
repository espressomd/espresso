#ifndef __MPI_CALLBACKS
#define __MPI_CALLBACKS

#include "communication.hpp"
#include <vector>
#include <functional>

class MpiCallbacks {
 public:
  typedef std::function<void (int)> Function;
  static int add(const Function &f) {
    const int id = m_callbacks.size();
    m_callbacks.push_back(f);
  }
  static void slave(int id, int par) {
    assert(id < m_callbacks.size());
    m_callbacks[id](par);    
  }
  static void call(int id) {
    mpi_call(MpiCallbacks::slave, id, 0);
  }
 private:
  static std::vector<Function> m_callbacks;
};

#endif
