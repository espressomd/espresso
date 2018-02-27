
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_tribend_helpers.hpp"
#include "utils/make_unique.hpp" //for creating a unique ptr to a bond class object
#include "bond/tBendingMethod.hpp" // for Bond::tBendingMethod::...



/****************
  IBM_Tribend_ResetParams
 *****************/

int IBM_Tribend_ResetParams(const int bond_type, const double kb)
{
  
  // Check if bond exists and is of correct type
  if ( bond_type >= n_bonded_ia ) {
    printf("bond does not exist while reading tribend checkpoint\n");
    return ES_ERROR;
  }

  auto ibm_tribend_bond = dynamic_cast<Bond::IbmTribend*>(bond_container.get_Bond(bond_type));
  if(ibm_tribend_bond){
    ibm_tribend_bond.ResetParams(kb);
  }
  else{
    printf("Dynamic Cast of Bond failed: interaction type does not match while reading tribend checkpoint!\n");
    return ES_ERROR; 
  };
  
  return ES_OK;
}

/***********
   IBM_Tribend_SetParams
************/

int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const int ind4, const tBendingMethod method, const double kb, const bool flat)
{

  
  // Distinguish bending methods
  if ( method == TriangleNormals )
  {
    double theta0;
  
    if ( !flat )
    {
      // Compute theta0
      //Particle p1, p2, p3, p4;
      auto p1 = get_particle_data(ind1);
      auto p2 = get_particle_data(ind2);
      auto p3 = get_particle_data(ind3);
      auto p4 = get_particle_data(ind4);
      
      //Get vectors of triangles
      double dx1[3], dx2[3], dx3[3];
      get_mi_vector(dx1, p1->r.p, p3->r.p);
      get_mi_vector(dx2, p2->r.p, p3->r.p);
      get_mi_vector(dx3, p4->r.p, p3->r.p);
      
      //Get normals on triangle; pointing outwards by definition of indices sequence
      double n1l[3], n2l[3];
      vector_product(dx1, dx2, n1l);
      vector_product(dx1, dx3, n2l);
      
      // Wolfgang here had a minus. It seems to work, so leave it in
      n2l[0] = -1*n2l[0];
      n2l[1] = -1*n2l[1];
      n2l[2] = -1*n2l[2];
      
      double n1[3], n2[3];
      unit_vector(n1l,n1);
      unit_vector(n2l,n2);
      
      
      //calculate theta by taking the acos of the scalar n1*n2
      double sc = scalar(n1,n2);
      if ( sc > 1.0) sc = 1.0;
      
      theta0 = acos(sc);
      double tmp[3];
      vector_product(n1,n2,tmp);
      
      const double desc = scalar(dx1,tmp);
      if ( desc < 0) theta0 = 2.0*PI-theta0;
      
    }
    else theta0 = 0;        // Flat


    // NOTE: This is the bare bending modulus used by the program.
    // If triangle pairs appear only once, the total bending force should get a factor 2
    // For the numerical model, a factor sqrt(3) should be added, see Gompper&Kroll J. Phys. 1996 and Krüger thesis
    // This is an approximation, it holds strictly only for a sphere

    bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::IbmTribend>
		     (kb, Bond::tBendingMethod::TriangleNormals, theta0, 3));
  }
  
  // Gompper
  if ( method == NodeNeighbors )
  {
    // Interpret ind2 as number of partners
    
    if ( ind1 != 5 && ind1 != 6) { printf("Gompper bending with %d partners seems strange. Are you sure?\n", ind2); return ES_ERROR; }
    
    bond_container.set_bond_by_type(bond_type, Utils::make_unique<Bond::IbmTribend>
    		     (kb, Bond::tBendingMethod::NodeNeighbors, 0.0, ind1));
  }
  
  return ES_OK;
  
}

#endif
