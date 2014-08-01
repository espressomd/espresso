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

#include "immersed-boundary/stretching_force_ibm.hpp"


#ifdef STRETCHING_FORCE_IMMERSED_BOUNDARY
#include "communication.hpp"

//default law neo-hookean use setmd stretching_force_ibm_law 1 to change to skalak
int stretching_force_law_ibm = 0;

int stretching_force_ibm_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka) {
	Particle part1, part2, part3;
	double templo[3], templpo[3],vecpro[3];
	double lo, lpo, sinpo, cospo;
	double a1, a2, a3, b1, b2, b3;
	double area2;
	
	if(bond_type < 0) {
    	return ES_ERROR;
    }

  	make_bond_type_exist(bond_type);
	
	//Get data (especially location) of three particles
	get_particle_data(ind1,&part1);
	get_particle_data(ind2,&part2);
	get_particle_data(ind3,&part3);
	
	//Calculate equilibrium lenghts and angle; Note the sequence of the points!
	//lo = lenght between 1 and 3
	//vecsub(part3.r.p,part1.r.p,templo);
	get_mi_vector(templo, part3.r.p, part1.r.p);
	lo = sqrt (sqrlen(templo));
	//lpo = lenght between 1 and 2
	//vecsub(part2.r.p,part1.r.p,templpo);
	get_mi_vector(templpo, part2.r.p, part1.r.p);
	lpo = sqrt (sqrlen(templpo));
	//cospo / sinpo angle functions between these vectors; calculated directly via the producs
	cospo = scalar(templo,templpo)/(lo*lpo);
	vector_product(templo, templpo, vecpro);
	sinpo = sqrt(sqrlen(vecpro))/(lo*lpo);
	
	//Use the values determined above for further constants of the stretch-force calculation
	area2 = lo * lpo * sinpo;
	a1 = -(lo*sinpo)/area2;
	a2 = - a1;
	a3 = 0.0;
	b1 = (lo*cospo - lpo)/area2;
	b2 = -(lo*cospo)/area2;
	b3 = lpo/area2;

	//Hand these values over to parameter structure
  	bonded_ia_params[bond_type].p.stretching_force_ibm.a1 = a1;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.a2 = a2;
	bonded_ia_params[bond_type].p.stretching_force_ibm.a3 = a3;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b1 = b1;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b2 = b2;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b3 = b3;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.lo = lo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.lpo = lpo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.sinpo = sinpo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.cospo = cospo;
	bonded_ia_params[bond_type].p.stretching_force_ibm.Area0 = 0.5*area2;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.maxdist = max;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.ks = ks;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.ka = ka;
  	
  	bonded_ia_params[bond_type].type = STRETCHING_FORCE_IBM_IA;
    bonded_ia_params[bond_type].num = 2;
  	
  	//Communicate this to whoever is interested
  	mpi_bcast_ia_params(bond_type, -1); 
	
	//Free particles's allocated memory 
	free_particle(&part1);
	free_particle(&part2);
	free_particle(&part3);	
	
	return ES_OK;
}

int stretching_force_ibm_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka) {
  
	double a1, a2, a3, b1, b2, b3;
	double area2;
	
	if(bond_type < 0) {
	  return ES_ERROR;
        }

  	make_bond_type_exist(bond_type);
	
	//Use the values determined above for further constants of the stretch-force calculation
	area2 = 2.0 * Area0;
	a1 = -(lo*sinpo)/area2;
	a2 = - a1;
	a3 = 0.0;
	b1 = (lo*cospo - lpo)/area2;
	b2 = -(lo*cospo)/area2;
	b3 = lpo/area2;

	//Hand these values over to parameter structure
	//bonded_ia_params[bond_type].p.stretching_force_ibm.ind1 = ind1;
	//bonded_ia_params[bond_type].p.stretching_force_ibm.ind2 = ind2;
	//bonded_ia_params[bond_type].p.stretching_force_ibm.ind3 = ind3;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.a1 = a1;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.a2 = a2;
	bonded_ia_params[bond_type].p.stretching_force_ibm.a3 = a3;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b1 = b1;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b2 = b2;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.b3 = b3;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.lo = lo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.lpo = lpo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.sinpo = sinpo;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.cospo = cospo;
	bonded_ia_params[bond_type].p.stretching_force_ibm.Area0 = 0.5*area2;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.maxdist = max;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.ks = ks;
  	bonded_ia_params[bond_type].p.stretching_force_ibm.ka = ka;
  	
  	bonded_ia_params[bond_type].type = STRETCHING_FORCE_IBM_IA;
    bonded_ia_params[bond_type].num = 2;
  	
  	//Communicate this to whoever is interested
  	mpi_bcast_ia_params(bond_type, -1); 

	
	return ES_OK;
}

#endif
