#include "cells.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "initialize.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "rotation.hpp"
#include "utils.hpp"

void local_rotate_system(double phi, double theta, double alpha) {
  // Culculate center of mass
  double com[3] = {0, 0, 0};
  int N = 0; // Num of particles

  for (auto const &p : local_cells.particles()) {
    for (int j = 0; j < 3; j++) {
      com[j] += p.r.p[j];
    }
    N++;
  }

  for (int j = 0; j < 3; j++)
    com[j] /= N;

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

  resort_particles = 1;
  announce_resort_particles();
}

/** Rotate all particle coordinates around an axis given by phi,theta through
 * the center of mass by an angle alpha */
void rotate_system(double phi, double theta, double alpha) {
  if (n_nodes != 1) {
    runtimeErrorMsg() << "Rotate_system only works on a single cpu core";
  }

  local_rotate_system(phi, theta, alpha);
}
