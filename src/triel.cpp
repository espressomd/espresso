#include "triel.h"


#ifdef TRIELASTIC
#include "communication.h"


int triel_set_params(int bond_type, int ind1, int ind2, int ind3, double max, double ks, double ka) {
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
  	bonded_ia_params[bond_type].p.triel.a1 = a1;
  	bonded_ia_params[bond_type].p.triel.a2 = a2;
	bonded_ia_params[bond_type].p.triel.a3 = a3;
  	bonded_ia_params[bond_type].p.triel.b1 = b1;
  	bonded_ia_params[bond_type].p.triel.b2 = b2;
  	bonded_ia_params[bond_type].p.triel.b3 = b3;
  	bonded_ia_params[bond_type].p.triel.lo = lo;
  	bonded_ia_params[bond_type].p.triel.lpo = lpo;
  	bonded_ia_params[bond_type].p.triel.sinpo = sinpo;
  	bonded_ia_params[bond_type].p.triel.cospo = cospo;
	bonded_ia_params[bond_type].p.triel.Area0 = 0.5*area2;
  	bonded_ia_params[bond_type].p.triel.maxdist = max;
  	bonded_ia_params[bond_type].p.triel.ks = ks;
  	bonded_ia_params[bond_type].p.triel.ka = ka;
  	
  	bonded_ia_params[bond_type].type = TRIEL_IA;
    bonded_ia_params[bond_type].num = 2;
  	
  	//Communicate this to whoever is interested
  	mpi_bcast_ia_params(bond_type, -1); 
	
	//Free particles's allocated memory 
	free_particle(&part1);
	free_particle(&part2);
	free_particle(&part3);	
	
	return ES_OK;
}

int triel_reset_params(int bond_type, double lo, double lpo, double cospo, double sinpo, double Area0, double max, double ks, double ka) {
  
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
	//bonded_ia_params[bond_type].p.triel.ind1 = ind1;
	//bonded_ia_params[bond_type].p.triel.ind2 = ind2;
	//bonded_ia_params[bond_type].p.triel.ind3 = ind3;
  	bonded_ia_params[bond_type].p.triel.a1 = a1;
  	bonded_ia_params[bond_type].p.triel.a2 = a2;
	bonded_ia_params[bond_type].p.triel.a3 = a3;
  	bonded_ia_params[bond_type].p.triel.b1 = b1;
  	bonded_ia_params[bond_type].p.triel.b2 = b2;
  	bonded_ia_params[bond_type].p.triel.b3 = b3;
  	bonded_ia_params[bond_type].p.triel.lo = lo;
  	bonded_ia_params[bond_type].p.triel.lpo = lpo;
  	bonded_ia_params[bond_type].p.triel.sinpo = sinpo;
  	bonded_ia_params[bond_type].p.triel.cospo = cospo;
	bonded_ia_params[bond_type].p.triel.Area0 = 0.5*area2;
  	bonded_ia_params[bond_type].p.triel.maxdist = max;
  	bonded_ia_params[bond_type].p.triel.ks = ks;
  	bonded_ia_params[bond_type].p.triel.ka = ka;
  	
  	bonded_ia_params[bond_type].type = TRIEL_IA;
    bonded_ia_params[bond_type].num = 2;
  	
  	//Communicate this to whoever is interested
  	mpi_bcast_ia_params(bond_type, -1); 

	
	return ES_OK;
}

#endif
