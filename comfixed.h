// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.

/** \file comfixed.h
 *  Routines to enable comfixed
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:sayar@mpip-mainz.mpg.de">Mehmet</a>
*/
#ifdef COMFIXED
MDINLINE void calc_comfixed()
{
  Particle *p;
  int i, np, c;
  Cell *cell;
  IA_parameters *ia_params;
  int t0,type_count0;
  int j;
	double fsum0[3];

  for (t0=0; t0<n_particle_types; t0++) {
    ia_params = get_ia_param(t0,t0);
    if(ia_params->COMFIXED_flag == 1) {
  	  type_count0=0;
      for(j = 0; j < 3; j++) {
			  fsum0[j]= 0.;
			}
      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
	 	      if(p[i].p.type==t0) {
			      type_count0 ++;
      	    for(j = 0; j < 3; j++) {
				      fsum0[j] += p[i].f.f[j];
			      }
		      }
        }
      }

      for (c = 0; c < local_cells.n; c++) {
        cell = local_cells.cell[c];
        p  = cell->part;
        np = cell->n;
        for(i = 0; i < np; i++) {
	 	      if(p[i].p.type==t0) {
      	    for(j = 0; j < 3; j++) {
				      p[i].f.f[j] -= fsum0[j]/type_count0;
			      }
		      }
        }
      }
    }
  }  
}
#endif
