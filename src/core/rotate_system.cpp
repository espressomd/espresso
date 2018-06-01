#include "cells.hpp"
#include "communication.hpp"
#include "particle_data.hpp"
#include "rotation.hpp"
#include "utils.hpp"
#include "initialize.hpp"

#include <boost/mpi/collectives.hpp>

namespace mpi = boost::mpi;

void local_rotate_system(double phi, double theta, double alpha) {
  // Culculate center of mass
  Vector3d local_com{};
  double local_mass = 0.0;

  for (auto const &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      local_com[j] += p.p.mass * p.r.p[j];
    }
    local_mass += p.p.mass;
  }

  auto const total_mass =
      mpi::all_reduce(comm_cart, local_mass, std::plus<double>());
  auto const com =
      mpi::all_reduce(comm_cart, local_com, std::plus<Vector3d>()) / total_mass;

  // Rotation axis in carthesian coordinates
  double axis[3];
  axis[0] = sin(theta) * cos(phi);
  axis[1] = sin(theta) * sin(phi);
  axis[2] = cos(theta);

  // Rotate particle coordinates
  for (auto &p : local_cells.particles()) {
    // Move the center of mass of the system to the origin
    for (int j = 0; j < 3; j++) {
      p.r.p[j] -= com[j];
    }
    // Rotate
    double res[3];
    vec_rotate(axis, alpha, p.r.p, res);
    // Write back result and shift back the center of mass
    for (int j = 0; j < 3; j++) {
      p.r.p[j] = com[j] + res[j];
    }
#ifdef ROTATION
    rotate_particle(&p, axis, alpha);
#endif
  }

  set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

void mpi_rotate_system_slave(int, int) {
  std::array<double, 3> params;
  mpi::broadcast(comm_cart, params, 0);

  local_rotate_system(params[0], params[1], params[2]);
}

void mpi_rotate_system(double phi, double theta, double alpha) {
  mpi_call(mpi_rotate_system_slave, 0, 0);

  std::array<double, 3> params{{phi, theta, alpha}};
  mpi::broadcast(comm_cart, params, 0);

  local_rotate_system(params[0], params[1], params[2]);
}

/** Rotate all particle coordinates around an axis given by phi,theta through
 * the center of mass by an angle alpha */
void rotate_system(double phi, double theta, double alpha) {
  mpi_rotate_system(phi, theta, alpha);
}
