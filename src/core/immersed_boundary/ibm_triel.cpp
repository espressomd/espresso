#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object


/****************
  IBM_Triel_ResetParams
 *****************/
int IBM_Triel_ResetParams(const int bond_type, const double k1, const double l0)
{
  
  // Check if bond exists and is of correct type
  if (bond_type >= n_bonded_ia) {
    printf("bond does not exist while reading triel checkpoint\n");
    return ES_ERROR;
  }

  auto ibm_triel_bond = dynamic_cast<Bond::IbmTriel*>(bond_container.get_Bond(bond_type));
  if(ibm_triel_bond){
    return ibm_triel_bond.ResetParams(k1, l0);
  }
  else{
    printf("Dynamic Cast of Bond failed: interaction type does not match while reading triel checkpoint!\n");
    return ES_ERROR;
  };
  
  return ES_OK;
}

/***********
   IBM_Triel_SetParams
************/
int IBM_Triel_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const double max, const tElasticLaw elasticLaw, const double k1, const double k2)
{

  //Get data (especially location) of three particles
  //Particle part1, part2, part3;
  auto part1 = get_particle_data(ind1);
  auto part2 = get_particle_data(ind2);
  auto part3 = get_particle_data(ind3);
  
  // Calculate equilibrium lenghts and angle; Note the sequence of the points!
  // lo = lenght between 1 and 3
  double templo[3];
  get_mi_vector(templo, part3->r.p, part1->r.p);
  const double l0 = sqrt (sqrlen(templo));
  // lpo = lenght between 1 and 2
  double templpo[3];
  get_mi_vector(templpo, part2->r.p, part1->r.p);
  const double lp0 = sqrt (sqrlen(templpo));

  // cospo / sinpo angle functions between these vectors; calculated directly via the products
  const double cosPhi0 = scalar(templo,templpo)/(l0*lp0);
  double vecpro[3];
  vector_product(templo, templpo, vecpro);
  const double sinPhi0 = sqrt(sqrlen(vecpro))/(l0*lp0);
  
  // Use the values determined above for further constants of the stretch-force calculation
  const double area2 = l0 * lp0 * sinPhi0;
  const double a1 = -(l0*sinPhi0)/area2;
  const double a2 = - a1;
  const double b1 = (l0*cosPhi0 - lp0)/area2;
  const double b2 = -(l0*cosPhi0)/area2;

  //create bond class
  // Always store two constants k1, k2, for NeoHookean only k1 is used
  bond_container.set_bond_by_type(bond_type,
				  Utils::make_unique<Bond::IbmTriel>(a1, a2, b1, b2, l0, lp0,
								     sinPhi0, cosPhi0, 0.5*area2,
								     max, elasticLaw, k1, k2));
  
  return ES_OK;
}
#endif
