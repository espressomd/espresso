/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file statistics_molecule.cpp

    see \ref statistics_molecule.hpp
*/
#include "statistics_molecule.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

/* new version for new topology structure */
int analyze_fold_molecules(float *coord, double shift[3])
{
  int m,p,i, tmp;
  int mol_size, ind;
  double cm_tmp, com[3];

  /* check molecule information */
  if ( n_molecules < 0 ) return ES_ERROR;

  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{059 analyze_fold_molecules: could not sort particle config, particle ids not consecutive?} ");
    return ES_ERROR;
  }

  /* loop molecules */
  for(m=0; m<n_molecules; m++) {
    mol_size = topology[m].part.n;
    if(mol_size > 0) {
      /* calc center of mass */
      calc_mol_center_of_mass(topology[m],com);
      /* fold coordinates */
      for(i=0; i<3; i++) {
	if ( PERIODIC(i) ) { 
	  tmp = (int)floor((com[i]+shift[i])*box_l_i[i]);
	  cm_tmp =0.0;
	  for(p=0; p<mol_size; p++) {
	    ind        = 3*topology[m].part.e[p] + i;
	    coord[ind] -= tmp*box_l[i];
	    coord[ind] += shift[i];
	    cm_tmp     += coord[ind];
	  }
	  cm_tmp /= (double)mol_size;
	  if(cm_tmp < -10e-6 || cm_tmp > box_l[i]+10e-6) {
	    char *errtxt = runtime_error(128 + ES_INTEGER_SPACE + 2*ES_DOUBLE_SPACE);
	    ERROR_SPRINTF(errtxt,"{060 analyze_fold_molecules: chain center of mass is out of range (coord %d: %.14f not in box_l %.14f)} ",
		    i,cm_tmp,box_l[i]);
	    return ES_ERROR;
	  }
	}
      }
    }
  }
  return ES_OK;
}



/* TODO: this function is not used anywhere. To be removed? */
double calc_mol_hydro_radius(Molecule mol) 
{
  int i, j, id1, id2;
  double rh=0.0, diff_vec[3];

  for(i=0; i<mol.part.n; i++) {
    id1 = mol.part.e[i];
    for(j=i+1; j<mol.part.n; j++) {
      id2 = mol.part.e[i];
      vecsub(partCfg[id1].r.p, partCfg[id2].r.p, diff_vec);
      rh += 1.0/sqrt(sqrlen(diff_vec));
    }
  }
  return 0.5*(mol.part.n*(mol.part.n-1))/rh;
}


/**Incorporates mass of each particle*/
void calc_mol_center_of_mass(Molecule mol, double com[3])
{
  int i,j,id;
  double M = 0.0;
  for(j=0; j<3; j++) com[j]=0.0;

  for(i=0; i<mol.part.n; i++) {
    id = mol.part.e[i];
    for(j=0; j<3; j++) com[j]+= partCfg[id].r.p[j]*PMASS(partCfg[id]);
    M += PMASS(partCfg[id]);
  }
    for(j=0; j<3; j++) com[j] /= M;
}

/* TODO: This function is not used anywhere. To be removed? */
/**Incorporates mass of each particle*/
double calc_mol_gyr_radius2(Molecule mol)
{
  int i, id;
  double rg=0.0, M=0.0, com[3], diff_vec[3];

  calc_mol_center_of_mass(mol, com);

  for(i=0; i<mol.part.n; i++) {
    id = mol.part.e[i];
    vecsub(partCfg[id].r.p, com, diff_vec);
    rg += sqrlen(diff_vec);
    M += PMASS(partCfg[id]);
  }

  return (rg/M);
}
