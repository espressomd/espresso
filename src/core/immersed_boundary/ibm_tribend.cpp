
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "immersed_boundary/ibm_tribend.hpp"

// DEBUG
/*double maxBendingForce;
double maxBendingDist;
double maxX;*/

/*************
   IBM_Tribend_CalcForce
Calculate the bending force and add it to the particles
 **************/

void IBM_Tribend_CalcForce(Particle *p1, const int numPartners, Particle **const partners, const Bonded_ia_parameters &iaparams)
{
    // move to separate function
    if ( numPartners != 3 ) { 
        throw std::runtime_error("IBM_Tribend expects 3 bond partners.");
    }   
    
    
    Particle *p2 = partners[0];
    Particle *p3 = partners[1];
    Particle *p4 = partners[2];

    assert(p1);
    assert(p2);
    assert(p3);
    assert(p4);
    
    // ************* This is Wolfgang's code **************
    // with some corrections by Achim
    
    //Get vectors making up the two triangles
    double dx1[3], dx2[3], dx3[3];
    get_mi_vector(dx1, p1->r.p, p3->r.p);
    get_mi_vector(dx2, p2->r.p, p3->r.p);
    get_mi_vector(dx3, p4->r.p, p3->r.p);
    
    //Get normals on triangle; pointing outwards by definition of indices sequence
    double n1[3], n2[3];
    vector_product(dx1, dx2, n1);
    vector_product(dx1, dx3, n2);
    
    // Wolfgang here had a minus. It seems to work, so leave it in
    n2[0]=-1*n2[0];
    n2[1]=-1*n2[1];
    n2[2]=-1*n2[2];
    
    
    // Get 2*area of triangles out of the magnitude of the resulting normals and make the latter unity
    const double Ai = sqrt(n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2]);
    n1[0] = n1[0]/Ai;
    n1[1]=n1[1]/Ai;
    n1[2]=n1[2]/Ai;
    
    const double Aj = sqrt(n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2]);
    n2[0] = n2[0]/Aj;
    n2[1]=n2[1]/Aj;
    n2[2]=n2[2]/Aj;
    
    
    //Get the prefactor for the force term
    double sc = scalar(n1,n2);
    if ( sc > 1.0) sc = 1.0;
    
    //Get theta as angle between normals
    double theta = acos(sc);
    
    double direc[3];
    vector_product(n1,n2,direc);
    const double desc = scalar(dx1,direc);
    
    if (desc < 0) theta = -1.0*theta;
    const double DTh = theta - iaparams.p.ibm_tribend.theta0;
   
    double Pre;
    // Correct version with linearized sin
    if ( theta > 0) Pre = 1.0*iaparams.p.ibm_tribend.kb * DTh;
    else Pre = -1.0*iaparams.p.ibm_tribend.kb * DTh;
    
    double v1l[3], v2l[3];
    for (int i = 0; i < 3; i++)
    {
      v1l[i] = n2[i]-sc*n1[i];
      v2l[i] = n1[i]-sc*n2[i];
    }
    
    double len = sqrt(sqrlen(v1l));
    double v1[3], v2[3];
    if (len>0) {
      for ( int i = 0;i <3; i++)  v1[i]=v1l[i]/len;
    }
    else
    {
       throw std::runtime_error("Case len<=0 not handled.");
    }
    
    // Achim normalizes both with the length of v1, Wolfgang uses v1 and v2
    // However, the length should be identical, so it does not matter

    if ( len > 0) for (int i = 0;i <3; i++) {
      v2[i]=v2l[i]/len;
    }
    else
    {
       throw std::runtime_error("Case len<=0 not handled.");
    }
    
    
    
    // Force for particle 1:
    double tmp[3], tmp2[3], term1[3], term2[3];
    get_mi_vector(tmp,p2->r.p,p3->r.p);
    get_mi_vector(tmp2, p3->r.p, p4->r.p);
    vector_product(tmp,v1, term1);
    vector_product(tmp2,v2, term2);
    for (int i = 0; i < 3; i++ )
      p1->f.f[i] += Pre*(term1[i]/Ai + term2[i]/Aj);
    
     // Force for particle 2:
    get_mi_vector(tmp,p3->r.p,p1->r.p);
    vector_product(tmp,v1, term1);
    for (int i = 0; i < 3; i++)
      p2->f.f[i] += Pre*(term1[i]/Ai);
    
    // Force for Particle 3:
    get_mi_vector(tmp,p1->r.p,p2->r.p);
    get_mi_vector(tmp2, p4->r.p, p1->r.p);
    vector_product(tmp,v1, term1);
    vector_product(tmp2,v2, term2);
    for (int i = 0; i < 3; i++)
      p3->f.f[i] += Pre*(term1[i]/Ai + term2[i]/Aj);
    
    // Force for Particle 4:
    get_mi_vector(tmp,p1->r.p,p3->r.p);
    vector_product(tmp,v2, term1);
    for (int i = 0; i < 3; i++) p4->f.f[i] += Pre*(term1[i]/Aj);
  }

/****************
  IBM_Tribend_ResetParams
 *****************/

int IBM_Tribend_ResetParams(const int bond_type, const double kb)
{
  
  // Check if bond exists and is of correct type
  if ( bond_type >= bonded_ia_params.size() ) { printf("bond does not exist while reading tribend checkpoint\n"); return ES_ERROR; }
  if ( bonded_ia_params[bond_type].type != BONDED_IA_IBM_TRIBEND ) { printf("interaction type does not match while reading tribend checkpoint!\n"); return ES_ERROR; }

  // Check if k is correct
  if ( fabs( bonded_ia_params[bond_type].p.ibm_tribend.kb - kb) > 1e-6 ) { printf("kb does not match while reading tribend checkpoint. It is %.12e and read was %.12e\n", bonded_ia_params[bond_type].p.ibm_tribend.kb, kb); return ES_ERROR; }
  
  //Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);
  
  return ES_OK;
}

/***********
   IBM_Tribend_SetParams
************/

int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const int ind4, const double kb, const bool flat)
{
  // Create bond
  make_bond_type_exist(bond_type);
  
  // General parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_TRIBEND;
  
  // Specific parameters
//  bonded_ia_params[bond_type].p.ibm_tribend.method = method;
  
  // Distinguish bending methods
//  if ( method == TriangleNormals )
  {
    double theta0;
  
    if ( !flat )
    {
      // Compute theta0
      auto p1 = get_particle_data(ind1);
      auto p2 = get_particle_data(ind2);
      auto p3 = get_particle_data(ind3);
      auto p4 = get_particle_data(ind4);
      
      //Get vectors of triangles
      double dx1[3], dx2[3], dx3[3];
      get_mi_vector(dx1, p1.r.p, p3.r.p);
      get_mi_vector(dx2, p2.r.p, p3.r.p);
      get_mi_vector(dx3, p4.r.p, p3.r.p);
      
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

    // Krüger always has three partners
    bonded_ia_params[bond_type].num = 3;
    bonded_ia_params[bond_type].p.ibm_tribend.theta0 = theta0;
    // NOTE: This is the bare bending modulus used by the program.
    // If triangle pairs appear only once, the total bending force should get a factor 2
    // For the numerical model, a factor sqrt(3) should be added, see Gompper&Kroll J. Phys. 1996 and Krüger thesis
    // This is an approximation, it holds strictly only for a sphere
    bonded_ia_params[bond_type].p.ibm_tribend.kb = kb;
  }

  // Broadcast and return
  mpi_bcast_ia_params(bond_type, -1);
  return ES_OK;
}

#endif
