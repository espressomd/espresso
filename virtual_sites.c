
// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.

#include "virtual_sites.h"
#include "pressure.h"

#ifdef VIRTUAL_SITES

// The following four functions are independent of the specif
// rules used to place virtual particles

void update_mol_vel_pos()
{
   //replace this by a better implementation later!

   // ORDER MATTERS! Update_mol_vel may rely on correct positions of virtual particcles
   update_mol_pos();
   update_mol_vel();
}

void update_mol_vel()
{
#ifndef VIRTUAL_SITES_NO_VELOCITY
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       if (ifParticleIsVirtual(&p[i]))
        update_mol_vel_particle(&p[i]);
    }
  }
#endif
}

void update_mol_pos()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
       if (ifParticleIsVirtual(&p[i]))
        update_mol_pos_particle(&p[i]);
    }
    //only for real particles
  }
}

int update_mol_pos_cfg(){
  int i;
  for(i=0; i<n_total_particles; i++) {
     if (partCfg[i].p.isVirtual==1)
      update_mol_pos_particle(&partCfg[i]);
  }
  return 1;
}

// Now, load the rules deciding, how to place virtual particles
// and transfer back forces to the real particles
#ifdef VIRTUAL_SITES_COM
 #include "virtual_sites_com.c"
#endif
#ifdef VIRTUAL_SITES_RELATIVE
 #include "virtual_sites_relative.c"
#endif

#endif 

