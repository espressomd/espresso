#ifndef __MPI_CALLBACKS
#define __MPI_CALLBACKS

#include <vector>
#include <functional>
#include <iostream>

void mpi_init(int *argc, char ***argv);
namespace boost {
namespace mpi {
class communicator;
}}

class MpiCallbacks {
 public:
  typedef std::function<void (int)> Function;
  
  static int add(const Function &f);  
  static void call(int id, int par);

  static boost::mpi::communicator mpi_comm;
  
 private:
  friend void mpi_init(int *argc, char ***argv);
  static void slave(int id, int par);
  static std::vector<Function> m_callbacks;
};

#endif
