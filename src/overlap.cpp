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
/** \file overlap.cpp
 *
 *  Implementation of \ref overlap.hpp
 */
#include "utils.hpp"
#include "overlap.hpp"
#include "interaction_data.hpp"
#include "communication.hpp"

#ifdef OVERLAPPED

int overlapped_bonded_set_params(int bond_type, int overlap_type,
				 char * filename) 
{
  int i, scan_success = 0, size;
  FILE* fp;

  if(bond_type < 0)
    return 1;
  
  make_bond_type_exist(bond_type);

  /* set types */
  bonded_ia_params[bond_type].type       = BONDED_IA_OVERLAPPED;
  bonded_ia_params[bond_type].p.overlap.type = overlap_type;

  /* set number of interaction partners */
  if(overlap_type == OVERLAP_BOND_LENGTH)   bonded_ia_params[bond_type].num = 1;
  if(overlap_type == OVERLAP_BOND_ANGLE) bonded_ia_params[bond_type].num = 2;
  if(overlap_type == OVERLAP_BOND_DIHEDRAL) bonded_ia_params[bond_type].num = 3;

  /* set max bondlength of between two beads, in Unit Angstrom */
  bonded_ia_params[bond_type].p.overlap.maxval = 6.0;

  /* copy filename */
  size = strlen(filename);
  bonded_ia_params[bond_type].p.overlap.filename = (char*)malloc((size+1)*sizeof(char));
  strcpy(bonded_ia_params[bond_type].p.overlap.filename,filename);

  fp = fopen( filename , "r");
  if ( !fp )
    return 2;
  
  /* Read in size of overlapps from file */
  scan_success = fscanf( fp , "%d ", &size);
  if ( scan_success < 1 ) { 
    fclose(fp);
    return 3;
  } 

  bonded_ia_params[bond_type].p.overlap.noverlaps = size;

  /* allocate overlapped funciton parameter arrays */
  bonded_ia_params[bond_type].p.overlap.para_a = (double*)malloc(size*sizeof(double));
  bonded_ia_params[bond_type].p.overlap.para_b = (double*)malloc(size*sizeof(double));
  bonded_ia_params[bond_type].p.overlap.para_c = (double*)malloc(size*sizeof(double));

   /* Read in the overlapped funciton parameter data */
  for (i=0; i<size; i++) {
  	scan_success = fscanf(fp, "%lg ", &bonded_ia_params[bond_type].p.overlap.para_a[i]);
  	if ( scan_success < 1 ) { 
    		fclose(fp);
    		return 3;
  	} 
  	scan_success = fscanf( fp, "%lg ", &bonded_ia_params[bond_type].p.overlap.para_b[i]);
  	if ( scan_success < 1 ) { 
    		fclose(fp);
    		return 3;
  	} 
	scan_success = fscanf( fp, "%lg ", &bonded_ia_params[bond_type].p.overlap.para_c[i]);
  	if ( scan_success < 1 ) { 
    		fclose(fp);
    		return 3;
  	} 
  }
  fclose(fp);

  mpi_bcast_ia_params(bond_type, -1); 

  return 0;
}

#endif

