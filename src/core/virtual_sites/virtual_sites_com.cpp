/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#include "virtual_sites_com.hpp"

#ifdef VIRTUAL_SITES_COM

#include "cells.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp"
#include "virtual_sites.hpp"

// forward declarations
Utils::Vector3d calc_mol_vel(Particle *p_com);
Utils::Vector3d calc_mol_pos(Particle *p_com);
void put_mol_force_on_parts(Particle *p_com);

void update_mol_vel_particle(Particle *p) {
  if (p->p.is_virtual) {
    p->m.v = calc_mol_vel(p);
  }
}

void update_mol_pos_particle(Particle *p) {
  if (p->p.is_virtual) {
    p->r.p = calc_mol_pos(p);
  }
}

void distribute_mol_force(const ParticleRange &particles) {
  for (auto &p : particles) {
    if (p.p.is_virtual) {
      if (p.f.f.norm2() != 0) {
        put_mol_force_on_parts(&p);
      }
    }
  }
}

Utils::Vector3d calc_mol_vel(Particle *p_com) {
  double M = 0;
  Particle *p;
  Utils::Vector3d v_com{};
  int const mol_id = p_com->p.mol_id;
  for (int i = 0; i < topology[mol_id].part.n; i++) {
    p = local_particles[topology[mol_id].part.e[i]];
    if (p->p.is_virtual)
      continue;
    v_com += p->p.mass * p->m.v;
    M += p->p.mass;
  }
  v_com /= M;
  return v_com;
}

/* this is a local version of center of mass, because ghosts don't have image
 * boxes*/
/* but p_com is a real particle */
Utils::Vector3d calc_mol_pos(Particle *p_com) {
  double M = 0;
  Particle *p;
  Utils::Vector3d r_com{};
  int const mol_id = p_com->p.mol_id;
  for (int i = 0; i < topology[mol_id].part.n; i++) {
    p = local_particles[topology[mol_id].part.e[i]];
    if (p->p.is_virtual)
      continue;
    r_com += p->p.mass * get_mi_vector(p->r.p, p_com->r.p, box_geo);
    M += p->p.mass;
  }
  r_com /= M;
  r_com += p_com->r.p;
  return r_com;
}

void put_mol_force_on_parts(Particle *p_com) {
  Particle *p;
  int const mol_id = p_com->p.mol_id;
  auto const force = p_com->f.f;
  p_com->f.f = {};
#ifdef MASS
  double M = 0;
  for (int i = 0; i < topology[mol_id].part.n; i++) {
    p = local_particles[topology[mol_id].part.e[i]];
    if (p->p.is_virtual)
      continue;
    M += p->p.mass;
  }
#else
  double M = topology[mol_id].part.n - 1;
#endif
  for (int i = 0; i < topology[mol_id].part.n; i++) {
    p = local_particles[topology[mol_id].part.e[i]];
    if (!p->p.is_virtual) {
      p->f.f += p->p.mass * force / M;
    }
  }
}

Particle *get_mol_com_particle(Particle *calling_p) {
  int const mol_id = calling_p->p.mol_id;

  if (mol_id < 0) {
    runtimeErrorMsg() << "Particle does not have a mol id! pnr= "
                      << calling_p->p.identity << "\n";
    return nullptr;
  }
  for (int i = 0; i < topology[mol_id].part.n; i++) {
    Particle *p = local_particles[topology[mol_id].part.e[i]];

    if (p == nullptr) {
      runtimeErrorMsg()
          << "Particle does not exist in get_mol_com_particle! id= "
          << topology[mol_id].part.e[i] << "\n";
      return nullptr;
    }

    if (p->p.is_virtual) {
      return p;
    }
  }

  runtimeErrorMsg() << "No COM found in get_mol_com_particle! pnr= "
                    << calling_p->p.identity << "\n";
  return nullptr;
}

#endif
