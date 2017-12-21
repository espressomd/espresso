
#include "config.hpp"

#ifdef IMMERSED_BOUNDARY

#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "grid.hpp"
#include "communication.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_tribend_helpers.hpp"

// DEBUG
/*double maxBendingForce;
double maxBendingDist;
double maxX;*/

// Internal function
void CalcForceGompper(Particle *p1, const int numPartners, Particle **const partners, const double kb);

/*************
   IBM_Tribend_CalcForce
Calculate the bending force and add it to the particles
 **************/

void IBM_Tribend_CalcForce(Particle *p1, const int numPartners, Particle **const partners, const Bonded_ia_parameters &iaparams)
{

  const tBendingMethod method = iaparams.p.ibm_tribend.method;
  if ( method == NodeNeighbors ) CalcForceGompper(p1, numPartners, partners, iaparams.p.ibm_tribend.kb);
  if ( method == TriangleNormals )
  {
    // move to separate function
    if ( numPartners != 3 ) { printf("Error. TriangleNormals bending with != 3 partners!\n"); exit(1); }
    Particle *p2 = partners[0];
    Particle *p3 = partners[1];
    Particle *p4 = partners[2];
    
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
    // Classical Wolfgang version
/*    if ( theta > 0) Pre = 1.0*iaparams.p.ibm_tribend.kb * sin(DTh);
      else Pre = -1.0*iaparams.p.ibm_tribend.kb * sin(DTh);*/
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
    if (len>0) for ( int i = 0;i <3; i++)  v1[i]=v1l[i]/len;
    
    // Achim normalizes both with the length of v1, Wolfgang uses v1 and v2
    // However, the length should be identical, so it does not matter
//    if ( method == Krueger )
//      len = sqrt(sqrlen(v2l));
    
    if ( len > 0) for (int i = 0;i <3; i++)  v2[i]=v2l[i]/len;
    
    
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
}

/****************
  IBM_Tribend_ResetParams
 *****************/

int IBM_Tribend_ResetParams(const int bond_type, const double kb)
{
  
  // Check if bond exists and is of correct type
  if ( bond_type >= n_bonded_ia ) { printf("bond does not exist while reading tribend checkpoint\n"); return ES_ERROR; }
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

int IBM_Tribend_SetParams(const int bond_type, const int ind1, const int ind2, const int ind3, const int ind4, const tBendingMethod method, const double kb, const bool flat)
{
  // Create bond
  make_bond_type_exist(bond_type);
  
  // General parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_TRIBEND;
  
  // Specific parameters
  bonded_ia_params[bond_type].p.ibm_tribend.method = method;
  
  // Distinguish bending methods
  if ( method == TriangleNormals )
  {
    double theta0;
  
    if ( !flat )
    {
      // Compute theta0
      Particle p1, p2, p3, p4;
      get_particle_data(ind1, &p1);
      get_particle_data(ind2, &p2);
      get_particle_data(ind3, &p3);
      get_particle_data(ind4, &p4);
      
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
  
  // Gompper
  if ( method == NodeNeighbors )
  {
    // Interpret ind2 as number of partners
    
    if ( ind1 != 5 && ind1 != 6) { printf("Gompper bending with %d partners seems strange. Are you sure?\n", ind2); return ES_ERROR; }
    
    bonded_ia_params[bond_type].num = ind1;
    
    // Only flat eq possible, but actually this is ignored in the computation anyway
    bonded_ia_params[bond_type].p.ibm_tribend.theta0 = 0;
    bonded_ia_params[bond_type].p.ibm_tribend.kb = kb;
  }
  
  // Broadcast and return
  mpi_bcast_ia_params(bond_type, -1);
  return ES_OK;
  
}

/**************
   CalcForceGompper
 **************/

void CalcForceGompper(Particle *xi, const int numNeighbors, Particle **const neighbors, const double kb)
{
//  if ( xi->r.p[0] > maxX ) maxX = xi->r.p[0];
  
  // DEBUG stuff
/*  for (int i=0; i < numNeighbors; i++)
  {
    Vector3D tmp;
    Subtr(tmp, xi, neighbors[i]);
    const double dist = Length(tmp);
    if ( dist > maxBendingDist ) maxBendingDist = dist;
  }*/
  
  // mainNumerator: Will be set to one of the terms in the numerator.
  Vector3D mainNumerator;
  mainNumerator.el[0] = mainNumerator.el[1] = mainNumerator.el[2] = 0;
  
  // Will be set to the denominator.
  double denominator = 0;
  
  // We want to differentiate with respect to the i'th node ("virtual" index numDerivatives-1) and with respect to the neighbours of the i'th node (indecies 0<=l<=numDerivatives-2).
  const int numDerivatives = numNeighbors + 1;
  
  // The derivatives of the numerator in the energy with respect to the i'th node and its neighbours.
  Matrix3D *derivNumerators = new Matrix3D[numDerivatives];
  
  // The derivatives of the denominator in the energy with respect to the i'th node and its neighbours.
  Vector3D *derivDenominators = new Vector3D[numDerivatives];
  
  for (int i=0; i < numDerivatives; i++)
  {
    for ( int k=0; k < 3; k++)
    {
      derivDenominators[i].el[k] = 0;
      for (int l=0; l < 3; l++)
        derivNumerators[i].el[k][l] = 0;
    }
  }
  
  // We need to calculate several sums over the neighbouring sites of i. This is done in the for-loop.
  for (int curNeighborIdx = 0; curNeighborIdx < numNeighbors; ++curNeighborIdx)
  {
    // Get the neighbours of xi.
    const Particle *xj = neighbors[curNeighborIdx];
    
    
    // Neighbor list with periodic mapping
    const Particle *xjM1 = (curNeighborIdx == 0) ? neighbors[numNeighbors-1] : neighbors[curNeighborIdx-1];
    const Particle *xjP1 = (curNeighborIdx == numNeighbors-1) ? neighbors[0] : neighbors[curNeighborIdx+1];
    
    Vector3D tmp;
    Subtr(tmp, xi, xj);
    const double dxijLen2 = LengthSqr(tmp);
    
    // Cosine of the angles between the neighbours.
    const double cosThetaM1 = CalcCosTheta(xi, xj, xjM1);
    const double cosThetaP1 = CalcCosTheta(xi, xj, xjP1);
    
    // The sum of the two cotangens.
    const double Tij = CalcCot(cosThetaM1) + CalcCot(cosThetaP1);
    
    Vector3D ijDifference;
    Subtr(ijDifference, xi, xj);
    
    // Update the two terms which are independent of the node xl (the derivatives are with respect to node xl).
    Multiply(tmp, ijDifference, Tij);
    AddTo( mainNumerator, tmp);
    denominator += dxijLen2 * Tij;
    
    
    // Now for the other terms, dependent on xl.
    // TODO: The only terms which are not 0 are jM1, j, jP1 and i
    for (int gradientNodeIdx = 0; gradientNodeIdx < numDerivatives; ++gradientNodeIdx)
    {
      
      const Particle *const xl = (gradientNodeIdx == numDerivatives-1) ? xi : neighbors[gradientNodeIdx];
      
      // Derivatives of the two cotangens.
      Vector3D TijDeriv;
      CalcCotDerivativeGompperAnalyt(TijDeriv, xi, xj, xjM1, xl->p.identity);
      Vector3D tmp;
      CalcCotDerivativeGompperAnalyt(tmp, xi, xj, xjP1, xl->p.identity);
      AddTo( TijDeriv, tmp);
      
      // Kronecker Deltas.
      const int ilDelta = (xi->p.identity == xl->p.identity) ? 1 : 0;
      const int jlDelta = (xj->p.identity == xl->p.identity) ? 1 : 0;
      
      // Update
      Matrix3D tmpM;
      DyadicProduct(tmpM, ijDifference, TijDeriv);
      AddTo( derivNumerators[gradientNodeIdx], tmpM );
      AddScalarTo( derivNumerators[gradientNodeIdx], Tij*(ilDelta-jlDelta) );
      
      Multiply(tmp, ijDifference, Tij*2. * (ilDelta-jlDelta) );
      AddTo( derivDenominators[gradientNodeIdx], tmp );
      Multiply(tmp, TijDeriv, dxijLen2);
      AddTo( derivDenominators[gradientNodeIdx], tmp );
    }
  }
  
  
  
  // Calculate the contribution to the force for each node.
  for (int gradientNodeIdx = 0; gradientNodeIdx < numDerivatives; ++gradientNodeIdx)
  {
    Particle *const xl = (gradientNodeIdx == numDerivatives-1) ? xi : neighbors[gradientNodeIdx];
    
    // -1: The force is minus the gradient of the energy.
    // Note: left vector-matrix product: mainNumerator * derivNumerators[gradientNodeIdx]
    Vector3D force, tmp;
    LeftVectorMatrix(force, mainNumerator, derivNumerators[gradientNodeIdx]);
    Multiply(force, force, 4.0 / denominator );
    Multiply(tmp, derivDenominators[gradientNodeIdx], - 2.0 * LengthSqr(mainNumerator) / (denominator*denominator));
    AddTo(force, tmp);
    
    Multiply(force, force, (-1.0) * (kb / 2.0));
    
//    if ( xl->p.identity == 0)
//      printf("  Adding to node = %d when treating node %d: f = %.12e %.12e %.12e\n", xl->p.identity, xi->p.identity, force.el[0], force.el[1], force.el[2]);
    xl->f.f[0] += force.el[0];
    xl->f.f[1] += force.el[1];
    xl->f.f[2] += force.el[2];
    
    // DEBUG stuff
/*    for (int i=0; i < numNeighbors; i++)
    {
      const double f = Length(force);
      if ( f > maxBendingForce ) maxBendingForce = f;
    }*/
    
  }
  
  // Clean up
  delete []derivNumerators;
  delete []derivDenominators;
    
}



#endif
