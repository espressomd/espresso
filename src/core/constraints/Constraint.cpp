#include <boost/mpi/collectives.hpp>

#include "Constraint.hpp"
#include "communication.hpp"
#include "energy_inline.hpp"
#include "errorhandling.hpp"
#include "forces_inline.hpp"
#include "interaction_data.hpp"

namespace Constraints {

Vector3d Constraint::total_force() const {
  Vector3d total_force;
  boost::mpi::all_reduce(comm_cart, m_local_force, total_force,
                         std::plus<Vector3d>());

  return total_force;
}

}
