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
    assert(f != nullptr);
    const int id = m_callbacks.size();
    std::cout << this_node << ": Adding callback id = " << id << std::endl;        
    m_callbacks.push_back(f);

    return id;
  }
  static void slave(int id, int par) {
    std::cout << this_node << ": MpiCallbacks::slave id = " << id << std::endl;        
    assert(id < m_callbacks.size());
    m_callbacks[id](par);
  }
  static void call(int id, int par) {
    std::cout << this_node << ": MpiCallbacks::call id = " << id << std::endl;        
    assert((id >= 0) && (id < m_callbacks.size()));
    mpi_call(MpiCallbacks::slave, id, par);
  }
 private:
  static std::vector<Function> m_callbacks;
};

#endif
