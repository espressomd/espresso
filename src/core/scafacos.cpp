/*
  Copyright (C) 2010,2011,2012,2013,2014,2015 The ESPResSo project
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

#include <vector>
#include <array>
#include <cassert>
#include <memory>
#include <string>
#include <algorithm>
#include <iostream>
#include <limits>

#include <fcs.h>
#include <mpi.h>

#include "cells.hpp"
#include "errorhandling.hpp"
#include "communication.hpp"
#include "initialize.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "tuning.hpp"
#include "interaction_data.hpp"

namespace Electrostatics {
namespace Scafacos {

struct ScafacosData {
  int update_particle_data();
  int update_particle_forces() const;
  /** Inputs */
  std::vector<double> charges, positions;
  /** Outputs */
  std::vector<double> forces, potentials;
};

#define handle_error(stmt) { const FCSResult res = stmt; if(res) runtimeError(fcs_result_get_message(res)); }

struct Scafacos : public ScafacosData {
  Scafacos(const std::string &_method, MPI_Comm comm, const std::string &parameters) : method(_method) {
    handle_error(fcs_init(&handle, method.c_str(), comm));

    int near_flag;
    fcs_get_near_field_delegation(handle, &near_flag);
    has_near = near_flag != 0;

    fcs_set_resort(handle, 0);
    
    parse_parameters(parameters);

    std::cout << "Scafacos() has_near = " << has_near << " parameters '" << parameters << "'" << std::endl;
  }
  ~Scafacos() {
    fcs_destroy(handle);
  }
  /** Set parameters common to all methods */
  void set_common_parameters(double *box_l, int *periodicity, int total_particles);
  double pair_force(double dist) const;
  void run();
  void tune();
  void parse_parameters(const std::string &s) {
    handle_error(fcs_parser(handle, s.c_str(), 0));
  }
  double r_cut() {
    if(has_near) {     
      fcs_float r_cut;
      
      fcs_get_r_cut(handle, &r_cut);      
      return r_cut;
    }
    else {
      return 0.0;
    }
  }
  
  void set_r_cut(double r_cut) {
    if(has_near) {
      fcs_set_r_cut(handle, r_cut);
    }
  }
 
  
  /** Handle from the library */
  FCS handle;
  /** Whether the method supports near field delegation */
  bool has_near;
  /** The scafacos method of this instance */
  const std::string method;
};

Scafacos *scafacos = 0;

void add_pair_force(Particle *p1, Particle *p2, double *d, double dist, double *force) {
  assert(scafacos);
  const double field = scafacos->pair_force(dist);
      
  for(int i=0; i<3; i++){
    p1->f.f[i] -= field * d[i] /dist * p2->p.q * p1->p.q * coulomb.bjerrum;
    p2->f.f[i] += field * d[i] /dist * p2->p.q * p1->p.q * coulomb.bjerrum;
  }
}

void add_long_range_force() {
  std::cout << "Electrostatics::Scafacos::add_long_range_force()" << std::endl;
  if(scafacos)
    scafacos->run();
  else
    runtimeError("Scafacos internal error.");
}


/** Determine runtime for a specific cutoff */
double time_r_cut(double r_cut) {
  double t;

  /** Set cutoff to time */
  scafacos->set_r_cut(r_cut);
  puts("set_r_cut()");
  /** Tune other parameters */
  scafacos->tune();
  puts("scafacos->tune()");
  
  on_coulomb_change();
  
  return time_force_calc(10);
}

/** Determine the optimal cutoff by bisection */
void tune_r_cut() {
  std::cout << "tune_r_cut()" << std::endl;
  const double tune_limit = 1e-3;
  /** scafacos p3m and ewald do not accept r_cut 0 for no good reason */
  double r_min = 1.0;
  double r_max = std::min(min_local_box_l, min_box_l/2.0) - skin;
  double t_min = 0;
  double t_max = std::numeric_limits<double>::max();

  /** Run bisection */
  while(std::fabs(r_min - r_max) > tune_limit) {
    std::cout << "Bisection loop with r_min " << r_min << " r_max " << r_max << std::endl;
    const double dr = 0.5*(r_max - r_min);
    std::cout << "t_mid" << std::endl;
    const double t_mid = time_r_cut(r_min + dr);
    std::cout << "t_min" << std::endl;
    t_min = time_r_cut(r_min);
    std::cout << "t_max" << std::endl;
    t_max = time_r_cut(r_max);

    if(t_min <= 0.0) {
      r_min += tune_limit;
      break;
    }

    if(t_max <= 0.0) {
      r_min -= tune_limit;
      break;
    }
    
    std::cout << "t_min " << t_min << " t_mid " << t_mid << " t_max " << t_max  << std::endl;
    
    if(t_mid > t_min) {
      r_max = r_min += dr;
    } else {
      r_min += dr;
    }    
  }  
}

void tune() {
  /** Check whether we have to do a bisection for the short range cutoff */
  /** Check if there is a user supplied cutoff */
  if((scafacos->has_near) && (scafacos->r_cut() <= 0.0)) {
    tune_r_cut();
  } else {
    scafacos->tune();
  }  
}

static void set_params_safe(const std::string &method, const std::string &params) {
  if(scafacos && (scafacos->method != method)) {
    delete scafacos;
    scafacos = 0;
  } else {
    scafacos = new Scafacos(method, comm_cart, params);
  }

  scafacos->parse_parameters(params);
  int per[3] = { PERIODIC(0) != 0, PERIODIC(1) != 0, PERIODIC(2) != 0 };
  scafacos->set_common_parameters(box_l, per, n_part);

  on_coulomb_change();
  
  tune();

  on_coulomb_change();
}

void set_parameters_slave(int n_method, int n_params) {
  std::string method;
  std::string params;

  method.resize(n_method);
  params.resize(n_params);
  
  /** This requires C++11, otherwise this is undefined because std::string was not required to have conitnuous memory before. */
  MPI_Bcast(&(*method.begin()), n_method, MPI_CHAR, 0, comm_cart);
  MPI_Bcast(&(*params.begin()), n_params, MPI_CHAR, 0, comm_cart);
    
  set_params_safe(method, params);
}

void set_parameters(const std::string &method, const std::string &params) {
  mpi_call(set_parameters_slave, method.size(), params.size());

  /** This requires C++11, otherwise this is undefined because std::string was not required to have conitnuous memory before. */
  /* const_cast is ok, this code runs only on rank 0 where the mpi call does not modify the buffer */
  MPI_Bcast(const_cast<char *>(&(*method.begin())), method.size(), MPI_CHAR, 0, comm_cart);
  MPI_Bcast(const_cast<char *>(&(*params.begin())), params.size(), MPI_CHAR, 0, comm_cart);

  set_params_safe(method, params);
}

double Scafacos::pair_force(double dist) const {
  if(has_near) {
    fcs_float field;
    fcs_compute_near_field(handle, dist, &field);
    return field;
  }

  return 0.0;
}

void Scafacos::run() {
  std::cout << "Electrostatics::Scafacos::Scafacos::run()" << std::endl;
  const int local_n_part = update_particle_data();
  std::cout << "local_n_part = " << local_n_part << std::endl;
  forces.resize(3*local_n_part);
  potentials.resize(local_n_part);
  std::cout << "resized." << std::endl;

  std::cout << "fcs_tune()...";
  handle_error(fcs_tune(handle, local_n_part, &(positions[0]), &(charges[0])));
  std::cout << "done.";
  
  //  handle_error(fcs_run(handle, local_n_part, &(positions[0]), &(charges[0]), &(forces[0]), &(potentials[0])));
  std::cout << "fcs_run()...";
  handle_error(fcs_run(handle, local_n_part, &(positions[0]), &(charges[0]), &(forces[0]), 0));
  std::cout << "done.";
  
  update_particle_forces();
}

void Scafacos::tune() {
  const int local_n_part = update_particle_data();
  
  handle_error(fcs_tune(handle, local_n_part, &(positions[0]), &(charges[0])));  
}

/** \brief Collect particle data in continous arrays as required by fcs */
int ScafacosData::update_particle_data() {
  std::cout << "ScafacosData::update_particle_data()" << std::endl;
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
int ScafacosData::update_particle_forces() const {
  std::cout << "forces.size() = " << forces.size() << std::endl;
  
 for(int c = 0; c < local_cells.n; c++) {
    Cell const * cell = local_cells.cell[c];
    Particle *p = cell->part;
    const int np = cell->n;
    
    for(int i = 0; i < np; i++) {
      p[i].f.f[0] += coulomb.bjerrum*p[i].p.q*forces[3*i + 0];
      p[i].f.f[1] += coulomb.bjerrum*p[i].p.q*forces[3*i + 1];
      p[i].f.f[2] += coulomb.bjerrum*p[i].p.q*forces[3*i + 2];
    }
 }
}

void Scafacos::set_common_parameters(double *box_l, int *periodicity, int total_particles) {
  double boxa[3] = { 0., 0., 0. }, boxb[3] = { 0., 0., 0. }, boxc[3] = { 0., 0., 0. }, off[3] = { 0., 0., 0. };
  boxa[0] = box_l[0];
  boxb[1] = box_l[1];
  boxc[2] = box_l[2];
  
  handle_error(fcs_set_common(handle, 0, boxa, boxb, boxc, off, periodicity, total_particles));

  //  fcs_print_parameters(handle);
}

double get_r_cut() {
  if(scafacos) {
    return scafacos->r_cut();
  }
  return 0.0;
}
}

}

