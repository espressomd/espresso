/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

#ifndef _MOL_CUT_H
#define _MOL_CUT_H

#include "virtual_sites.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"

#ifdef MOL_CUT

#define CUTOFF_CHECK(cond) ((ia_params->mol_cut_type!=0)||(cond))

///
int molcut_set_params(int part_type_a, int part_type_b,int mol_cut_type,double mol_cut_cutoff);

inline int checkIfParticlesInteractViaMolCut(Particle *p1, Particle *p2,IA_parameters *data){
   if (data->mol_cut_type==0){
      return 1;
   }
   else if (p1->p.mol_id == p2->p.mol_id) {
      return 1;
   }
   else{
      double com_dist=get_mol_dist(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

inline int checkIfParticlesInteractViaMolCut_partcfg(Particle *p1, Particle *p2,IA_parameters *data){
   if (data->mol_cut_type==0){
      return 1;
   }
   else if (p1->p.mol_id == p2->p.mol_id) {
      return 1;
   }
   else{
      double com_dist=get_mol_dist_partcfg(p1,p2);
      if (com_dist < data->mol_cut_cutoff) return 1;
   }
   return 0;
}

#else
#define CUTOFF_CHECK(cond) (cond)
#endif


#endif
