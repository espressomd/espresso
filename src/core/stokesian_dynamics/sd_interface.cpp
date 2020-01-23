/*
 * Copyright (C) 2010-2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.hpp"

#ifdef STOKESIAN_DYNAMICS

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include <utils/Vector.hpp>


#include <boost/config.hpp>

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp" // still needed?
#include "particle_data.hpp"
#include "ParticleRange.hpp"
#include "thermostat.hpp"

#if defined(BLAS) && defined(LAPACK)
#include "sd_cpu.hpp"
#endif

#ifdef CUDA
#include "sd_gpu.hpp"
#endif

#include "sd_interface.hpp"

namespace {
double sd_viscosity = -1.0;

enum { CPU, GPU, INVALID } device = INVALID;

std::unordered_map<int, double> radius_dict;

double sd_kT = 0.0;

std::size_t sd_seed = 0UL;

int sd_flags = 0;


/* Buffer that holds local particle data, and all particles on the master node 
   used for sending particle data to master node */
std::vector<SD_particle_data> parts_buffer{};

/* Buffer that holds the (translational and angular) velocities of the local
   particles on each node, used for returning results */
std::vector<double> v_sd{};


int debug_cnt = 0; 

} // namespace






/* packs selected properties of all local particles into a buffer */
std::vector<SD_particle_data> &sd_gather_local_particles(ParticleRange const &parts) {
  std::size_t n_part = parts.size();
  parts_buffer.resize(n_part);
  std::size_t i = 0;
  
  for (auto const &p: parts) {
    parts_buffer[i].id = p.p.is_virtual ? -1 : p.p.identity;
    parts_buffer[i].type = p.p.type;
    
    parts_buffer[i].pos[0] = p.r.p[0];
    parts_buffer[i].pos[1] = p.r.p[1];
    parts_buffer[i].pos[2] = p.r.p[2];
    
    // Velocity is not needed for SD, thus omit while gathering
    
    // TODO: rotation
    parts_buffer[i].quat[0] = 1.0;
    parts_buffer[i].quat[1] = 0.0;
    parts_buffer[i].quat[2] = 0.0;
    parts_buffer[i].quat[3] = 0.0;
    
    parts_buffer[i].ext_force[0] = p.f.f[0];
    parts_buffer[i].ext_force[1] = p.f.f[1];
    parts_buffer[i].ext_force[2] = p.f.f[2];
    
    parts_buffer[i].ext_torque[0] = p.f.torque[0];
    parts_buffer[i].ext_torque[1] = p.f.torque[1];
    parts_buffer[i].ext_torque[2] = p.f.torque[2];
    
    //TODO: radius?
    parts_buffer[i].r = radius_dict[p.p.type];
    //parts_buffer[i].r = 1;
    
    i++;
  }
  
  return parts_buffer;
}

std::vector<double> &get_sd_local_v_buffer() {
  v_sd.resize(local_cells.particles().size() * 6);
  return v_sd;
}

void sd_update_locally(ParticleRange const &parts) {
  std::size_t n_part = parts.size();
  std::size_t i = 0;
  
  for (auto &p: parts) {
    // Copy velocities
    p.m.v[0] = v_sd[6 * i];
    p.m.v[1] = v_sd[6 * i + 1];
    p.m.v[2] = v_sd[6 * i + 2];
    
    p.m.omega[0] = v_sd[6 * i + 3];
    p.m.omega[1] = v_sd[6 * i + 4];
    p.m.omega[2] = v_sd[6 * i + 5];
    
    // Do the integration
    // Question: is the time step the same on all nodes? Should be ...
    p.r.p[0] += p.m.v[0] * time_step;
    p.r.p[1] += p.m.v[1] * time_step;
    p.r.p[2] += p.m.v[2] * time_step;    
    
    // TODO: rotation
    
    i++;
  }
}


void set_sd_viscosity(double eta) { sd_viscosity = eta; }

double get_sd_viscosity() { return sd_viscosity; }

void set_sd_device(std::string const &dev) {
  if (dev == "cpu") {
    device = CPU;
  } else if (dev == "gpu") {
    device = GPU;
  } else {
    device = INVALID;
  }
}

std::string get_sd_device() {
  switch (device) {
  case CPU:
    return "cpu";
  case GPU:
    return "gpu";
  default:
    return "invalid";
  }
}

void set_sd_radius_dict(std::unordered_map<int, double> const &x) {
  radius_dict = x;
}

std::unordered_map<int, double> get_sd_radius_dict() { return radius_dict; }

void set_sd_kT(double kT) { sd_kT = kT; }

double get_sd_kT() { return sd_kT; }

void set_sd_seed(std::size_t seed) { sd_seed = seed; }

std::size_t get_sd_seed() { return sd_seed; }

void set_sd_flags(int flg) { sd_flags = flg; }

int get_sd_flags() { return sd_flags; }

void propagate_vel_pos_sd() {
  if (thermo_switch & THERMO_SD) {
    if (BOOST_UNLIKELY(sd_viscosity < 0.0)) {
      runtimeErrorMsg() << "sd_viscosity has an invalid value: " +
                               std::to_string(sd_viscosity);
      return;
    }

    if (BOOST_UNLIKELY(sd_kT < 0.0)) {
      runtimeErrorMsg() << "sd_kT has an invalid value: " +
                               std::to_string(sd_viscosity);
      return;
    }
    
    
    printf("Doing it 1\n");
    
    
    mpi_sd_gather_particles();
    
    std::size_t n_part = local_cells.particles().size(); // parts_buffer.size();
    
    printf("Doing it 2, n_part = %lu \n",n_part);
    
    parts_buffer.resize(n_part);
    
    static std::vector<int> id{};  // wanted these static
    static std::vector<double> x_host{};
    static std::vector<double> f_host{};
    static std::vector<double> a_host{};
    
    id.resize(n_part);
    x_host.resize(6 * n_part);
    f_host.resize(6 * n_part);
    a_host.resize(6 * n_part);
    
    printf("Doing it 3\n");
    
    std::size_t i = 0;
    for (auto const &p : parts_buffer) {  
      id[i] = parts_buffer[i].id; // p.p.is_virtual ? -1 : p.p.identity;
      
      x_host[6 * i + 0] = p.pos[0]; // p.r.p[0];
      x_host[6 * i + 1] = p.pos[1]; // p.r.p[1];
      x_host[6 * i + 2] = p.pos[2]; // p.r.p[2];
      // TODO: Rotation
      x_host[6 * i + 3] = 0; 
      x_host[6 * i + 4] = 0;
      x_host[6 * i + 5] = 0;
      
      printf("Particle %lu @ (%f,%f,%f)\n",i,p.pos[0],p.pos[1],p.pos[2]);
      
      f_host[6 * i + 0] = p.ext_force[0]; // p.f.f[0];
      f_host[6 * i + 1] = p.ext_force[1]; // p.f.f[1];
      f_host[6 * i + 2] = p.ext_force[2]; // p.f.f[2];
      f_host[6 * i + 3] = p.ext_torque[0]; // p.f.torque[0];
      f_host[6 * i + 4] = p.ext_torque[1]; // p.f.torque[1];
      f_host[6 * i + 5] = p.ext_torque[2]; // p.f.torque[2];
      
      double radius = p.r; // radius_dict[p.p.type];
      
      //if (debug_cnt<5)
        // position
        //printf("Particle %lu: x=(%f,%f,%f)\n",i,x_host[6 * i + 0],x_host[6 * i + 1],x_host[6 * i + 2]);
        // force
        //printf("Particle %lu: f=(%f,%f,%f)\n",i,f_host[6 * i + 0],f_host[6 * i + 1],f_host[6 * i + 2]);
        // id
        //printf("Particle %lu: id=%d\n",i,id[i]);
      
      
      if (BOOST_UNLIKELY(radius < 0.0)) {
        runtimeErrorMsg() << "particle of type " + std::to_string(p.type) + 
                                 " has an invalid radius: " +
                                 std::to_string(radius);
        return;
      }

      
      
      a_host[i] = radius;

      ++i;
    } 
    
    printf("Doing it 4\n");
    
    std::size_t offset = std::round(sim_time / time_step);
    //std::vector<double> v_sd;
    switch (device) {

#if defined(BLAS) && defined(LAPACK)
    case CPU:
      v_sd = sd_cpu(x_host, f_host, a_host, n_part, sd_viscosity,
                    sd_kT / time_step, offset, sd_seed, sd_flags);
      break;
#endif


#ifdef CUDA
    case GPU:
      v_sd = sd_gpu(x_host, f_host, a_host, n_part, sd_viscosity,
                    sd_kT / time_step, offset, sd_seed, sd_flags);
      break;
#endif

    default:
      runtimeErrorMsg()
          << "Invalid device for Stokesian dynamics. Available devices:"
#if defined(BLAS) && defined(LAPACK)
             " cpu"
#endif
#ifdef CUDA
             " gpu"
#endif
          ;
      return;
    }
    
    printf("Doing it 5\n");
    
    /*
    
    // put results into buffer
    for (int i = 0; i < n_part; ++i) {
      // skip virtual particles
      if (parts_buffer[i].id == -1) {
        continue;
      }

      SD_particle_data &p = parts_buffer[i];
      
      // velocities
      p.vel[0] = v_sd[6 * i + 0];
      p.vel[1] = v_sd[6 * i + 1];
      p.vel[2] = v_sd[6 * i + 2];
      p.omega[0] = v_sd[6 * i + 3];
      p.omega[1] = v_sd[6 * i + 4];
      p.omega[2] = v_sd[6 * i + 5];
      // integrated position
      p.pos[0] = x_host[6 * i + 0] + time_step * v_sd[6 * i + 0];
      p.pos[1] = x_host[6 * i + 1] + time_step * v_sd[6 * i + 1];
      p.pos[2] = x_host[6 * i + 2] + time_step * v_sd[6 * i + 2];
      // TODO: integrated rotation
      p.quat[0] = 1;
      p.quat[1] = 0;
      p.quat[2] = 0;
      p.quat[3] = 0;
    }    
    
    */
    
    // send results over to nodes
    mpi_sd_scatter_results();
    
    printf("Doing it 6\n");
    
    // update local buffer ... does it work? 
    sd_update_locally(local_cells.particles());  
    
    
    printf("Doing it 7\n");
    
    
    
    /*
    // prepare buffers to be sent via MPI: counting and resizing
    
    for (std::pair<int,SD_node_buffer> &element: node_buffers) {
      element.second.counter = 0;
    }
    
    for (auto const &p : parts_buffer) {
      // skip virtual particles
      if (p.id == -1) {
        continue;
      }
      
      ++node_buffers[p.on_node].counter;
    }
    
    for (std::pair<int,SD_node_buffer> &element: node_buffers) {
      // hopefully, the size doesn't change often
      element.second.buffer.resize(node.second.counter);
      element.second.counter = 0;
    }
    
    
    // sort particles by their nodes and put data into appropriate buffer
    for (int i = 0; i < n_part; ++i) {
      // skip virtual particles
      if (parts_buffer[i].id == -1) {
        continue;
      }
      
      // put results into buffer
      SD_node_buffer &current = node_buffers[p.on_node];
      auto &p = current.buffer[current.counter];
      
      p.node = parts_buffer[i].node;
      p.id   = parts_buffer[i].id;
      // velocities
      p.vel[0] = v_sd[6 * i + 0];
      p.vel[1] = v_sd[6 * i + 1];
      p.vel[2] = v_sd[6 * i + 2];
      p.omega[0] = v_sd[6 * i + 3];
      p.omega[1] = v_sd[6 * i + 4];
      p.omega[2] = v_sd[6 * i + 5];
      // integrated position
      p.pos[0] = x_host[6 * i + 0] + time_step * v_sd[6 * i + 0];
      p.pos[1] = x_host[6 * i + 1] + time_step * v_sd[6 * i + 1];
      p.pos[2] = x_host[6 * i + 2] + time_step * v_sd[6 * i + 2];
      // TODO: integrated rotation
      p.quat[0] = 1;
      p.quat[1] = 0;
      p.quat[2] = 0;
      p.quat[3] = 0;
      
      ++current.counter;
    }
    
    mpi_sd_scatter_results();
    
    */
    
    // Ende meine Version

  

  }
}

#endif // STOKESIAN_DYNAMICS
