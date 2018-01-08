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
#include "virtual_sites_com.hpp"

#ifdef VIRTUAL_SITES_COM

#include "virtual_sites.hpp"
#include "cells.hpp"
#include "topology.hpp"
#include "forces.hpp"
#include "partCfg_global.hpp"
#include "integrate.hpp"

// forward declarations
void calc_mol_vel(Particle *p_com,double v_com[3]);
void calc_mol_pos(Particle *p_com,double r_com[3]);
void put_mol_force_on_parts(Particle *p_com);

void update_mol_vel_particle(Particle *p){
   int j;
   double v_com[3];
   if (p->p.isVirtual) {
      calc_mol_vel(p,v_com);
      for (j=0;j<3;j++){
         p->m.v[j] = v_com[j];
      }
   }
}


void update_mol_pos_particle(Particle *p){
   int j;
   double r_com[3];
   if (p->p.isVirtual) {
      calc_mol_pos(p,r_com);
      for (j=0;j<3;j++){
         p->r.p[j] = r_com[j];
      }
   }
}

void distribute_mol_force() {
  for (auto &p : local_cells.particles()) {
    if (p.p.isVirtual) {
      if (sqrlen(p.f.f) != 0) {
        put_mol_force_on_parts(&p);
      }
    }
  }
}

void calc_mol_vel(Particle *p_com,double v_com[3]){
   int i,j,mol_id;
   double M=0;
   Particle *p;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
   for (i=0;i<3;i++){
      v_com[i]=0.0;
   }
   mol_id=p_com->p.mol_id;
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];
      #ifdef VIRTUAL_SITES_DEBUG
      if (p==nullptr){
          runtimeErrorMsg() <<"Particle does not exist in calc_mol_vel! id= " << topology[mol_id].part.e[i] << "\n";
         return;
      }
      #endif
      if (p->p.isVirtual) continue;
      for (j=0;j<3;j++){
         v_com[j] += (*p).p.mass*p->m.v[j];
      }
      M+=(*p).p.mass;
#ifdef VIRTUAL_SITES_DEBUG
      count++;
#endif
   }
   for (j=0;j<3;j++){
      v_com[j] /= M;
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
       runtimeErrorMsg() <<"There is more than one COM in calc_mol_vel! mol_id= " << mol_id << "\n";
      return;
   }
#endif
}

/* this is a local version of center of mass, because ghosts don't have image boxes*/
/* but p_com is a real particle */
void calc_mol_pos(Particle *p_com,double r_com[3]){
   int i,j,mol_id;
   double M=0;
   double vec12[3];
   Particle *p;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
   for (i=0;i<3;i++){
      r_com[i]=0.0;
   }
   mol_id=p_com->p.mol_id;
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];
      #ifdef VIRTUAL_SITES_DEBUG
      if (p==nullptr){
          runtimeErrorMsg() <<"Particle does not exist in calc_mol_pos! id= " << topology[mol_id].part.e[i] << "\n";
         return;
      }
      #endif
      if (p->p.isVirtual) continue;
      get_mi_vector(vec12,p->r.p, p_com->r.p);
      for (j=0;j<3;j++){
          r_com[j] += (*p).p.mass*vec12[j];
      }
      M+=(*p).p.mass;
#ifdef VIRTUAL_SITES_DEBUG
      count++;
#endif
   }
   for (j=0;j<3;j++){
      r_com[j] /= M;
      r_com[j] += p_com->r.p[j];
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
       runtimeErrorMsg() <<"There is more than one COM in calc_mol_pos! mol_id= " << mol_id << "\n";
      return;
   }
#endif
}

void put_mol_force_on_parts(Particle *p_com){
   int i,j,mol_id;
   Particle *p;
   double force[3],M;
#ifdef VIRTUAL_SITES_DEBUG
   int count=0;
#endif
  mol_id=p_com->p.mol_id;
   for (i=0;i<3;i++){
      force[i]=p_com->f.f[i];
      p_com->f.f[i]=0.0;
   }
#ifdef MASS
   M=0;
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];
      #ifdef VIRTUAL_SITES_DEBUG
      if (p==nullptr){
          runtimeErrorMsg() <<"Particle does not exist in put_mol_force_on_parts! id= " << topology[mol_id].part.e[i] << "\n";
         return;
      }
      #endif
       if (p->p.isVirtual) continue;
      M+=(*p).p.mass;
   }
#else
   M=topology[mol_id].part.n-1;
#endif
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];
      #ifdef VIRTUAL_SITES_DEBUG
      if (p==nullptr){
          runtimeErrorMsg() <<"Particle does not exist in put_mol_force_on_parts! id= " << topology[mol_id].part.e[i] << "\n";
         return;
      }
      #endif
      if (!p->p.isVirtual) {
         for (j=0;j<3;j++){
            p->f.f[j]+=(*p).p.mass*force[j]/M;
         }
#ifdef VIRTUAL_SITES_DEBUG
         count++;
#endif
      }
   }
#ifdef VIRTUAL_SITES_DEBUG
   if (count!=topology[mol_id].part.n-1){
       runtimeErrorMsg() <<"There is more than one COM input_mol_force_on_parts! mol_id= " << mol_id << "\n";
      return;
   }
#endif
}

Particle *get_mol_com_particle(Particle *calling_p){
   int mol_id;
   int i;
   Particle *p;

   mol_id=calling_p->p.mol_id;

   if (mol_id < 0) {
     runtimeErrorMsg() <<"Particle does not have a mol id! pnr= " << calling_p->p.identity << "\n";
     return nullptr;
   }
   for (i=0;i<topology[mol_id].part.n;i++){
      p=local_particles[topology[mol_id].part.e[i]];

      if (p==nullptr){
          runtimeErrorMsg() <<"Particle does not exist in put_mol_force_on_parts! id= " << topology[mol_id].part.e[i] << "\n";
         return nullptr;
      }

      if (p->p.isVirtual) {
          return p;
       }
   }

   runtimeErrorMsg() <<"No com found in get_mol_com_particleParticle does not exist in put_mol_force_on_parts! pnr= " << calling_p->p.identity << "\n";
   return nullptr;

   return calling_p;
}


#endif
