// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file statistics_molecule.c 

    see \ref statistics_molecule.h    

    <b>Responsible:</b>
    <a href="mailto:hanjo@mpip-mainz.mpg.de">Hanjo</a>
*/
#include "statistics_molecule.h"

/* new version for new topology structure */
int analyze_fold_chains(float *coord)
{
  int m,p,i, tmp; 
  int mol_size, ind;
  double cm_tmp, com[3];

  /* check molecule information */
  if ( n_molecules < 0 ) return (TCL_ERROR);

  /* loop molecules */
  for(m=0; m<n_molecules; m++) {
    mol_size = molecules[m].part.n;
    if(mol_size > 0) {
      /* calc center of mass */
      /*      com[0] = 0.0; com[1] = 0.0; com[2] = 0.0;
      for(p=0; p<mol_size; p++) {
	ind = 3*molecules[m].part.e[p];
	for(i=0;i<3;i++) com[i] += coord[ind+i];
      }		   
      for(i=0;i<3;i++) com[i] /= (double)mol_size; */
      calc_center_of_mass(molecules[m],com);
      /* fold coordinates */
      for(i=0;i<3;i++) {
	if ( PERIODIC(i) ) { 
	  tmp = (int)floor(com[i]*box_l_i[i]);
	  cm_tmp =0.0;
	  for(p=0; p<mol_size; p++) {
	    ind        = 3*molecules[m].part.e[p] + i;
	    coord[ind] -= tmp*box_l[i];
	    cm_tmp     += coord[ind];
	  }
	  cm_tmp /= (double)mol_size;
	  if(cm_tmp < -10e-6 || cm_tmp > box_l[i]+10e-6) {
	    fprintf(stderr,"\n: analyse_fold_chains: chain center of mass is out of range (coord %d: %.14f not in box_l %.14f), exiting\n",i,cm_tmp,box_l[i]);
	    errexit();
	  }
	}
      }
    }
  }
  return (TCL_OK); 
}

void calc_center_of_mass(Molecule mol, double com[3])
{
  int i,j,id;
  for(j=0;j<3;j++) com[j]=0.0;

  for(i=0;i<mol.part.n;i++) {
    id = mol.part.e[i];
    for(j=0;j<3;j++) com[j]+= partCfg[id].r.p[j];
  }

  for(j=0;j<3;j++) com[j] /= mol.part.n;

}
