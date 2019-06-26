#ifndef ESPRESSO_NODEGRID_HPP
#define ESPRESSO_NODEGRID_HPP

#include <boost/mpi/communicator.hpp>

#include <utils/Vector.hpp>

#include <boost/mpi/exception.hpp>



class NodeGrid {
public:
  NodeGrid(boost::mpi::communicator const &comm, Utils::Vector3i const &dims)
      : m_cart_comm(
            Utils::Mpi::cart_create(comm, dims, /* periodicity */ {1, 1, 1})) {}
  NodeGrid(boost::mpi::communicator const &comm)
      : NodeGrid(comm, Utils::Mpi::dims_create<3>(comm)) {}

private:
  boost::mpi::communicator m_cart_comm;

  /** The number of nodes in each spatial dimension. */
  Utils::Vector3i m_node_grid;
  /** position of node in node grid */
  Utils::Vector3i m_node_pos;
  /** the six nearest neighbors of a node in the node grid. */
  Utils::Vector<int, 6> m_node_neighbors;
};

#endif // ESPRESSO_NODEGRID_HPP
