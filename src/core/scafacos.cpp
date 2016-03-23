/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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

#include "scafacos.hpp"

#if defined(SCAFACOS) and defined(ELECTROSTATICS)

#include <vector>

#include <cassert>
#include <memory>
#include <algorithm>
#include <iostream>
#include <limits>

#include "cells.hpp"
#include "errorhandling.hpp"
#include "communication.hpp"
#include "initialize.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "tuning.hpp"
#include "interaction_data.hpp"
#include "electrostatics/scafacos/Scafacos.hpp"

/** This file contains the c-like interface for Scafacos */

namespace Electrostatics {
namespace Scafacos {

/** Get available scafacos methods */
std::list<std::string> available_methods() {
  return Scafacos::available_methods();
}

/** Encapsulation for the particle data needed by scafacos */
struct ScafacosData {
  int update_particle_data();
  void update_particle_forces() const;
  /** Inputs */
  std::vector<double> charges, positions;
  /** Outputs */
  std::vector<double> forces;
};


/** \brief Collect particle data in continous arrays as required by fcs */
int ScafacosData::update_particle_data() {
  positions.clear();
  charges.clear();

  for(int c = 0; c < local_cells.n; c++) {
    Cell const * cell = local_cells.cell[c];
    const Particle *p = cell->part;
    const int np = cell->n;
    
    for(int i = 0; i < np; i++) {
      positions.push_back(p[i].r.p[0]);
      positions.push_back(p[i].r.p[1]);
      positions.push_back(p[i].r.p[2]);
      charges.push_back(p[i].p.q);
    }    
  }
  
  return charges.size();
}

/** \brief Write forces back to particles */

void ScafacosData::update_particle_forces() const {
  int it = 0;
  
 for(int c = 0; c < local_cells.n; c++) {
    Cell const * cell = local_cells.cell[c];
    Particle *p = cell->part;
    const int np = cell->n;
    
    for(int i = 0; i < np; i++) {
      p[i].f.f[0] += coulomb.prefactor*p[i].p.q*forces[it++];
      p[i].f.f[1] += coulomb.prefactor*p[i].p.q*forces[it++];
      p[i].f.f[2] += coulomb.prefactor*p[i].p.q*forces[it++];
    }
 }

 /** Check that the particle number did not change */
 assert(it == forces.size());
}

static Scafacos *scafacos = 0;
static ScafacosData particles;

void add_pair_force(Particle *p1, Particle *p2, double *d, double dist, double *force) {
  assert(scafacos);
  const double field = scafacos->pair_force(dist);
  const double fak = p2->p.q * p1->p.q * field * coulomb.prefactor / dist;
  
  for(int i=0; i<3; i++){
    p1->f.f[i] -= fak * d[i];
    p2->f.f[i] += fak * d[i];
  }
}

void add_long_range_force() {
  particles.update_particle_data();
  
  if(scafacos)
    scafacos->run(particles.charges, particles.positions, particles.forces);
  else
    runtimeError("Scafacos internal error.");
  
  particles.update_particle_forces();
}

/** Determine runtime for a specific cutoff */
double time_r_cut(double r_cut) {
  double t;

  /** Set cutoff to time */
  scafacos->set_r_cut(r_cut);

  /** Tune other parameters */
  scafacos->tune(particles.charges, particles.positions);
  
  on_coulomb_change();
  
  return time_force_calc(10);
}

/** Determine the optimal cutoff by bisection */
void tune_r_cut() {
  const double tune_limit = 1e-3;
  /** scafacos p3m and ewald do not accept r_cut 0 for no good reason */
  double r_min = 1.0;
  double r_max = std::min(min_local_box_l, min_box_l/2.0) - skin;
  double t_min = 0;
  double t_max = std::numeric_limits<double>::max();

  /** Run bisection */
  while(std::fabs(r_min - r_max) > tune_limit) {
    const double dr = 0.5*(r_max - r_min);
    const double t_mid = time_r_cut(r_min + dr);
    t_min = time_r_cut(r_min);
    t_max = time_r_cut(r_max);

    if(t_min <= 0.0) {
      r_min += tune_limit;
      break;
    }

    if(t_max <= 0.0) {
      r_min -= tune_limit;
      break;
    }

    if(t_mid > t_min) {
      r_max = r_min += dr;
    } else {
      r_min += dr;
    }    
  }  
}

void tune() {
  particles.update_particle_data();
 
  /** Check whether we have to do a bisection for the short range cutoff */
  /** Check if there is a user supplied cutoff */
  if((scafacos->has_near) && (scafacos->r_cut() <= 0.0)) {
    tune_r_cut();
  } else {
    scafacos->tune(particles.charges, particles.positions);
  }
}

static void set_params_safe(const std::string &method, const std::string &params) {
  if(scafacos && (scafacos->method != method)) {
    delete scafacos;
    scafacos = 0;
  }
  
  scafacos = new Scafacos(method, comm_cart, params);

  scafacos->parse_parameters(params);
  int per[3] = { PERIODIC(0) != 0, PERIODIC(1) != 0, PERIODIC(2) != 0 };
  scafacos->set_common_parameters(box_l, per, n_part);

  on_coulomb_change();
  
  tune();

  on_coulomb_change();
}

/** Bend result from scafacos back to original format */
std::string get_parameters() {
  if(!scafacos) {
    return std::string();
  }

  std::string p = scafacos->get_parameters();
  
  std::replace(p.begin(), p.end(), ',', ' ');
  
  return p;
}

double get_r_cut() {
  if(scafacos) {
    return scafacos->r_cut();
  }
  return 0.0;  
}

void set_parameters(const std::string &method, const std::string &params) {
  mpi_call(mpi_scafacos_set_parameters_slave, method.size(), params.size());

  /** This requires C++11, otherwise this is undefined because std::string was not required to have conitnuous memory before. */
  /* const_cast is ok, this code runs only on rank 0 where the mpi call does not modify the buffer */
  MPI_Bcast(const_cast<char *>(&(*method.begin())), method.size(), MPI_CHAR, 0, comm_cart);
  MPI_Bcast(const_cast<char *>(&(*params.begin())), params.size(), MPI_CHAR, 0, comm_cart);

  set_params_safe(method, params);
}

}
}

#endif /* SCAFACOS */

void mpi_scafacos_set_parameters_slave(int n_method, int n_params) {
  #if defined(SCAFACOS) and defined(ELECTROSTATICS)
  using namespace Electrostatics::Scafacos;
  std::string method;
  std::string params;

  method.resize(n_method);
  params.resize(n_params);
  
  /** This requires C++11, otherwise this is undefined because std::string was not required to have conitnuous memory before. */
  MPI_Bcast(&(*method.begin()), n_method, MPI_CHAR, 0, comm_cart);
  MPI_Bcast(&(*params.begin()), n_params, MPI_CHAR, 0, comm_cart);
    
  set_params_safe(method, params);
  #endif /* SCAFACOS */
}
