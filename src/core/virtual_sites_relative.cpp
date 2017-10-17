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

#include "virtual_sites_relative.hpp"
#include "rotation.hpp"

#ifdef VIRTUAL_SITES_RELATIVE

using std::ostringstream;

// This is the "relative" implementation for virtual sites.
// Virtual particles are placed relative to the position of a real particle

// Obtain the real particle from which a virtual particle derives it's position
// Note: for now, we use the mol_di property of Particle
Particle* vs_relative_get_real_particle(Particle* p)
{
 return local_particles[p->p.vs_relative_to_particle_id];
}

// Update the pos of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_pos_particle(Particle *p)
{
 // First obtain the real particle responsible for this virtual particle:
 // Find the 1st real particle in the topology for the virtual particle's mol_id
 Particle *p_real = vs_relative_get_real_particle(p);
 // Check, if a real particle was found
 if (!p_real)
 {
     ostringstream msg;
     msg <<"virtual_sites_relative.cpp - update_mol_pos_particle(): No real particle associated with virtual site.\n";
     runtimeError(msg);
   return;
 }
 
 // Calculate the quaternion defining the orientation of the vecotr connectinhg
 // the virtual site and the real particle
 // This is obtained, by multiplying the quaternion representing the director
 // of the real particle with the quaternion of the virtual particle, which 
 // specifies the relative orientation.
 double q[4];
 multiply_quaternions(p_real->r.quat,p->p.vs_relative_rel_orientation,q);
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
  new_pos[i] =p_real->r.p[i] +director[i]/l*p->p.vs_relative_distance;
  // Handle the case that one of the particles had gone over the periodic
  // boundary and its coordinate has been folded
  if (PERIODIC(i)) 
  {
    tmp =p->r.p[i] -new_pos[i];
    //printf("%f\n",tmp);
    if (tmp > box_l[i]/2.) {
     //printf("greater than box_l/2 %f\n",tmp);
     p->r.p[i] =new_pos[i] + box_l[i];
    }
    else if (tmp < -box_l[i]/2.) {
     //printf("smaller than box_l/2 %f\n",tmp);
     p->r.p[i] =new_pos[i] - box_l[i];
    }
    else p->r.p[i] =new_pos[i];
   }
   else p->r.p[i] =new_pos[i];
 }
}

// Update the vel of the given virtual particle as defined by the real 
// particles in the same molecule
void update_mol_vel_particle(Particle *p)
{
 // NOT TESTED!
 
 // First obtain the real particle responsible for this virtual particle:
 Particle *p_real = vs_relative_get_real_particle(p);
 // Check, if a real particle was found
 if (!p_real)
 {
     ostringstream msg;
     msg <<"virtual_sites_relative.cpp - update_mol_pos_particle(): No real particle associated with virtual site.\n";
     runtimeError(msg);
   return;
 }
 
 // Calculate the quaternion defining the orientation of the vecotr connectinhg
 // the virtual site and the real particle
 // This is obtained, by multiplying the quaternion representing the director
 // of the real particle with the quaternion of the virtual particle, which 
 // specifies the relative orientation.
 double q[4];
 multiply_quaternions(p_real->r.quat,p->p.vs_relative_rel_orientation,q);
 // Calculate the director resulting from the quaternions
 double director[3];
 convert_quat_to_quatu(q,director);
 // normalize
 double l =sqrt(sqrlen(director));
 // Division comes in the loop below

 // Get omega of real particle in space-fixed frame
 double omega_space_frame[3];
 convert_omega_body_to_space(p_real,omega_space_frame);
 // Obtain velocity from v=v_real particle + omega_real_particle \times director
 vector_product(omega_space_frame,director,p->m.v);

 int i;
 // Add prefactors and add velocity of real particle
 for (i=0;i<3;i++)
 {
  // Scale the velocity by the distance of virtual particle from the real particle
  // Also, espresso stores not velocity but velocity * time_step
  p->m.v[i] *= time_step * p->p.vs_relative_distance/l;
  // Add velocity of real particle
  p->m.v[i] += p_real->m.v[i];
 }
}

// Distribute forces that have accumulated on virtual particles to the
// associated real particles
void distribute_mol_force() {
  // Iterate over all the particles in the local cells
  for (auto &p : local_cells.particles()) {
    // We only care about virtual particles
    if (ifParticleIsVirtual(&p)) {
      update_mol_pos_particle(&p);

      // First obtain the real particle responsible for this virtual particle:
      Particle *p_real = vs_relative_get_real_particle(&p);

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
      int j;
      //       printf("Particle %d gets torque from %f %f %f of particle
      //       %d\n",p_real->p.identity, p.f.f[0], p.f.f[1],p.f.f[2],
      //       p.p.identity);
      for (j = 0; j < 3; j++) {
        p_real->f.torque[j] += tmp[j];
        //	 printf("%f ",tmp[j]);
        p_real->f.f[j] += p.f.f[j];
      }
    }
  }
}

// Setup the virtual_sites_relative properties of a particle so that the given virtaul particle will follow the given real particle
int vs_relate_to(int part_num, int relate_to)
{
    // Get the data for the particle we act on and the one we wnat to relate
    // it to.
    auto p_current = get_particle_data(part_num);
    auto p_relate_to = get_particle_data(relate_to);
    if (!p_current || !p_relate_to) {
        ostringstream msg;
        msg <<"Could not retrieve particle data for the given id";
        runtimeError(msg);
      return ES_ERROR;
    }
    
    // get teh distance between the particles
    double d[3];
    get_mi_vector(d, p_current->r.p,p_relate_to->r.p);
    
    
    
    // Check, if the distance between virtual and non-virtual particles is larger htan minimum global cutoff
    // If so, warn user
    double l=sqrt(sqrlen(d));
    if (l>min_global_cut) {
        ostringstream msg;
        msg << "Warning: The distance between virtual and non-virtual particle (" << l << ") is\nlarger than the minimum global cutoff (" << min_global_cut << "). This may lead to incorrect simulations\nunder certain conditions. Set the \"System()\" class property \"min_global_cut\" to\nincrease the minimum cutoff.\n";
        runtimeWarning(msg);
      return ES_ERROR;
    }

    // Now, calculate the quaternions which specify the angle between 
    // the director of the particel we relate to and the vector
    // (paritlce_we_relate_to - this_particle)
    double quat[4];
    // The vs_relative implemnation later obtains the direcotr by multiplying
    // the quaternions representing the orientation of the real particle
    // with those in the virtual particle. The re quulting quaternion is then
    // converted to a director.
    // Whe have quat_(real particle) *quat_(virtual particle) 
    // = quat_(obtained from desired director)
    // Resolving this for the quat_(virtaul particle)

    //Normalize desired director
    int i;
    
    // If the distance between real & virtual particle is 0
    // we just set the relative orientation to 1 0 0 0, as it is irrelevant but
    // needs to be a valid quaternion
    if (l!=0)
    {
      for (i=0;i<3;i++)
        d[i]/=l;

      // Obtain quaternions from desired director
      double quat_director[4];
      convert_quatu_to_quat(d, quat_director);
  
      // Define quat as described above:
      double x=0;
      for (i=0;i<4;i++)
       x+=p_relate_to->r.quat[i]*p_relate_to->r.quat[i];
  
      quat[0]=0;
      for (i=0;i<4;i++)
       quat[0] +=p_relate_to->r.quat[i]*quat_director[i];
      
      quat[1] =-quat_director[0] *p_relate_to->r.quat[1] 
         +quat_director[1] *p_relate_to->r.quat[0]
         +quat_director[2] *p_relate_to->r.quat[3]
         -quat_director[3] *p_relate_to->r.quat[2];
      quat[2] =p_relate_to->r.quat[1] *quat_director[3] 
        + p_relate_to->r.quat[0] *quat_director[2] 
        - p_relate_to->r.quat[3] *quat_director[1] 
        - p_relate_to->r.quat[2] * quat_director[0];
      quat[3] =quat_director[3] *p_relate_to->r.quat[0]
        - p_relate_to->r.quat[3] *quat_director[0] 
        + p_relate_to->r.quat[2] * quat_director[1] 
        - p_relate_to->r.quat[1] *quat_director[2];
      for (i=0;i<4;i++)
       quat[i]/=x;
     
     
     // Verify result
     double qtemp[4];
     multiply_quaternions(p_relate_to->r.quat,quat,qtemp);
     for (i=0;i<4;i++)
       if (fabs(qtemp[i]-quat_director[i])>1E-9)
         fprintf(stderr, "vs_relate_to: component %d: %f instead of %f\n",
  	       i, qtemp[i], quat_director[i]);
    }
    else
    { 
     quat[0]=1;
     quat[1]=quat[2]=quat[3]=0;
    }

    // Set the particle id of the particle we want to relate to, the distnace
    // and the relative orientation
    if (set_particle_vs_relative(part_num, relate_to, l, quat) == ES_ERROR) {
      ostringstream msg;
      msg << "setting the vs_relative attributes failed";
      runtimeError(msg);
      return ES_ERROR;
    }
    set_particle_virtual(part_num,1);
   
   return ES_OK;
}


// Rigid body conribution to scalar pressure and stress tensor
void vs_relative_pressure_and_stress_tensor(double* pressure, double* stress_tensor)
{
  // Division by 3 volume is somewhere else. (pressure.cpp after all presure calculations)
  // Iterate over all the particles in the local cells

  for (auto &p : local_cells.particles()) {
    if (!ifParticleIsVirtual(&p))
      continue;

    update_mol_pos_particle(&p);

    // First obtain the real particle responsible for this virtual particle:
    Particle *p_real = vs_relative_get_real_particle(&p);

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

