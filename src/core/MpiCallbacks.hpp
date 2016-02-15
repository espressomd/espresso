#ifndef __MPI_CALLBACKS
#define __MPI_CALLBACKS

#include "communication.hpp"
#include <vector>
#include <functional>
#include <iostream>

class MpiCallbacks {
 public:
  typedef std::function<void (int)> Function;
  static int add(const Function &f) {
    const int id = m_callbacks.size();
    std::cout << this_node << ": Adding callback id = " << id << std::endl;        
    m_callbacks.push_back(f);

    return id;
  }
  static void slave(int id, int par) {
    assert(id < m_callbacks.size());
    m_callbacks[id](par);
  }
  static void call(int id, int par) {
    assert(id < m_callbacks.size());
    mpi_call(MpiCallbacks::slave, id, par);
    m_callbacks[id](par);
  }
 private:
  static std::vector<Function> m_callbacks;
};

#endif
