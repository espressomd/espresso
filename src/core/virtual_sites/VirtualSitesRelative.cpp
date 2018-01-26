/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include "config.hpp"
#include "VirtualSitesRelative.hpp"
#include "rotation.hpp"
#include "errorhandling.hpp" 

#ifdef VIRTUAL_SITES_RELATIVE

void VirtualSitesRelative::update(bool recalc_positions) const {

for (auto& p : local_cells.particles()) {
  if (!p.p.isVirtual) continue;

  if (recalc_positions)
    update_pos(p);

  if (have_velocity())
    update_vel(p);

}

}





// This is the "relative" implementation for virtual sites.
// Virtual particles are placed relative to the position of a real particle

// Obtain the real particle from which a virtual particle derives it's position
// Note: for now, we use the mol_di property of Particle

// Update the pos of the given virtual particle as defined by the real 
// particles in the same molecule
void VirtualSitesRelative::update_pos(Particle& p) const
{
 // First obtain the real particle responsible for this virtual particle:
 // Find the 1st real particle in the topology for the virtual particle's mol_id
 const Particle *p_real = local_particles[p.p.vs_relative_to_particle_id];
 // Check, if a real particle was found
 if (!p_real)
 {
     runtimeErrorMsg() <<"virtual_sites_relative.cpp - update_mol_pos_particle(): No real particle associated with virtual site.\n";
   return;
 }
 
 // Calculate the quaternion defining the orientation of the vecotr connectinhg
 // the virtual site and the real particle
 // This is obtained, by multiplying the quaternion representing the director
 // of the real particle with the quaternion of the virtual particle, which 
 // specifies the relative orientation.
 double q[4];
 multiply_quaternions(p_real->r.quat,p.p.vs_relative_rel_orientation,q);
 // Calculate the director resulting from the quaternions
 double director[3];
 convert_quat_to_quatu(q,director);
 // normalize
 double l =sqrt(sqrlen(director));
 // Division comes in the loop below
 
 // Calculate the new position of the virtual sites from
 // position of real particle + director 
 int i;
 double new_pos[3];
 double tmp;
 for (i=0;i<3;i++)
 {
  new_pos[i] =p_real->r.p[i] +director[i]/l*p.p.vs_relative_distance;
  double old=p.r.p[i];
  // Handle the case that one of the particles had gone over the periodic
  // boundary and its coordinate has been folded
  if (PERIODIC(i)) 
  {
    tmp =p.r.p[i] -new_pos[i];
    if (tmp > box_l[i]/2.) {
     p.r.p[i] =new_pos[i] + box_l[i];
    }
    else if (tmp < -box_l[i]/2.) {
     p.r.p[i] =new_pos[i] - box_l[i];
    }
    else p.r.p[i] =new_pos[i];
   }
   else p.r.p[i] =new_pos[i];
  // Has the vs moved by more than a skin
  if (fabs(old- p.r.p[i]) >skin) {
    runtimeErrorMsg() << "Virtual site "<<p.p.identity<< " has moved by more than the skin."<<old<<"->" <<p.r.p[i]; 
  }
  
 }
}

// Update the vel of the given virtual particle as defined by the real 
// particles in the same molecule
void VirtualSitesRelative::update_vel(Particle& p) const
{
 // First obtain the real particle responsible for this virtual particle:
 Particle *p_real = local_particles[p.p.vs_relative_to_particle_id];
 // Check, if a real particle was found
 if (!p_real)
 {
     runtimeErrorMsg() <<"virtual_sites_relative.cpp - update_mol_pos_particle(): No real particle associated with virtual site.\n";
   return;
 }
 
 double d[3];
 get_mi_vector(d,p.r.p,p_real->r.p);

 // Get omega of real particle in space-fixed frame
 double omega_space_frame[3];
 convert_omega_body_to_space(p_real,omega_space_frame);
 // Obtain velocity from v=v_real particle + omega_real_particle \times director
 vector_product(omega_space_frame,d,p.m.v);

 int i;
 // Add prefactors and add velocity of real particle
 for (i=0;i<3;i++)
 {
  // Scale the velocity by the distance of virtual particle from the real particle
  // Also, espresso stores not velocity but velocity * time_step
  p.m.v[i] *= time_step;
  // Add velocity of real particle
  p.m.v[i] += p_real->m.v[i];
 }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void VirtualSitesRelative::back_transfer_forces_and_torques() const {
  // Iterate over all the particles in the local cells
  for (auto &p : local_cells.particles()) {
    // We only care about virtual particles
    if (p.p.isVirtual) {
      // First obtain the real particle responsible for this virtual particle:
      Particle *p_real = local_particles[p.p.vs_relative_to_particle_id];

      // Get distance vector pointing from real to virtual particle, respecting
      // periodic boundary i
      // conditions
      double d[3];
      get_mi_vector(d, p.r.p, p_real->r.p);

      // The rules for transfering forces are:
      // F_realParticle +=F_virtualParticle
      // T_realParticle +=f_realParticle \times
      // (r_virtualParticle-r_realParticle)

      // Calculate torque to be added on real particle
      double tmp[3];
      vector_product(d, p.f.f, tmp);

      // Add forces and torques
      for (int j = 0; j < 3; j++) {
        p_real->f.torque[j] += tmp[j];
        p_real->f.f[j] += p.f.f[j];
      }
    }
  }
}


// Setup the virtual_sites_relative properties of a particle so that the given virtaul particle will follow the given real particle



// Rigid body conribution to scalar pressure and stress tensor
void VirtualSitesRelative::pressure_and_stress_tensor_contribution(double* pressure, double* stress_tensor) const {
  // Division by 3 volume is somewhere else. (pressure.cpp after all presure calculations)
  // Iterate over all the particles in the local cells

  for (auto &p : local_cells.particles()) {
    if (!p.p.isVirtual)
      continue;

    update_pos(p);

    // First obtain the real particle responsible for this virtual particle:
    const Particle *p_real = local_particles[p.p.vs_relative_to_particle_id];

    // Get distance vector pointing from real to virtual particle, respecting
    // periodic boundary i
    // conditions
    double d[3];
    get_mi_vector(d, p_real->r.p, p.r.p);

    // Stress tensor conribution
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        stress_tensor[k * 3 + l] += p.f.f[k] * d[l];

    // Pressure = 1/3 trace of stress tensor
    // but the 1/3 is applied somewhere else.
    *pressure += (p.f.f[0] * d[0] + p.f.f[1] * d[1] + p.f.f[2] * d[2]);
  }

  *pressure /= 0.5 * time_step * time_step;
  for (int i = 0; i < 9; i++)
    stress_tensor[i] /= 0.5 * time_step * time_step;

}

#endif

