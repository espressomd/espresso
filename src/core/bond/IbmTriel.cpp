#include "IbmTriel.hpp"
#include "grid.hpp" //for get_mit_vector
#include "interaction_data.hpp" // for NeoHookean

//force
int Bond::IbmTriel::calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, 
						    double force[3], double force2[3]) const {

  
  // Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
  // l = length between 1 and 3
  // get_mi_vector is an Espresso function which considers PBC
  double vec2[3] = {0., 0., 0.};
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  const double l = sqrt (sqrlen(vec2));
  
  // lp = lenght between 1 and 2
  double vec1[3] = {0., 0., 0.};
  get_mi_vector(vec1, p2->r.p, p1->r.p);
  const double lp = sqrt (sqrlen(vec1));

  // angles between these vectors; calculated directly via the products
  const double cosPhi = scalar(vec1,vec2)/(lp*l);
  double vecpro[3] = {0., 0., 0.};
  vector_product(vec1, vec2, vecpro);
  const double sinPhi = sqrt(sqrlen(vecpro))/(l*lp);
  
  // Check for sanity
  if(    (lp-m_lp0 > m_maxdist)
     ||  (l-m_l0 > m_maxdist))
  {
    return 1;
  }
  
  // Variables in the reference state
  const double l0 = m_l0;
  const double lp0 = m_lp0;
  const double cosPhi0 = m_cosPhi0;
  const double sinPhi0 = m_sinPhi0;
  const double a1 = m_a1;
  const double a2 = m_a2;
  const double b1 = m_b1;
  const double b2 = m_b2;
  const double A0 = m_area0;
  
  
  // Displacement gradient tensor D: Eq. (C.9) Krüger thesis
  const double Dxx = lp/lp0;
  const double Dxy = ((l/l0*cosPhi) - (lp/lp0*cosPhi0)) / sinPhi0;
  const double Dyx = 0.0;
  const double Dyy = l/l0*sinPhi / sinPhi0;
  
  // Tensor G: (C.12)
  const double Gxx = Utils::sqr(Dxx)+Utils::sqr(Dyx);
  const double Gxy = Dxx*Dxy + Dyx*Dyy;
  const double Gyx = Dxx*Dxy + Dyy*Dyx;   // = Gxy because of symmetry
  const double Gyy = Utils::sqr(Dxy) + Utils::sqr(Dyy);
  
  // Strain invariants, C.11 and C.12
  const double i1 = (Gxx + Gyy) - 2;
  const double i2 = ((Gxx * Gyy) - (Gxy * Gyx)) - 1;
  
  // Derivatives of energy density E used in chain rule below: Eq. (C.14)
  double dEdI1;
  double dEdI2;
  if ( m_elasticLaw == NeoHookean )
  {
    // Neo-Hookean
    dEdI1 = m_k1/6.0;
    dEdI2 = (-1)*m_k1/(6.0*(i2+1.0)*(i2+1.0));
  }
  else
  {
    //Skalak
    dEdI1 = m_k1*(i1+1)/6.0;
    dEdI2 = (-1)*m_k1/6.0 + m_k2*i2/6.0;
  }
  
  // ******** Achim's version *****************
  
  // Derivatives of Is (C.15)
  const double dI1dGxx = 1;
  const double dI1dGxy = 0;
  const double dI1dGyx = 0;
  const double dI1dGyy = 1;
  
  const double dI2dGxx = Gyy;
  const double dI2dGxy = -Gyx; // Note: Krueger has a factor 2 here, because he uses the symmetry of the G-matrix.
  const double dI2dGyx = -Gxy; //			But we don't use it. So, Krueger is missing the yx term, whereas we have it.
  const double dI2dGyy = Gxx;
  
  // Derivatives of G (C.16)
  const double dGxxdV1x = 2*a1*Dxx;
  const double dGxxdV1y = 0;
  const double dGxxdV2x = 2*a2*Dxx;
  const double dGxxdV2y = 0;
  
  const double dGxydV1x = a1*Dxy + b1*Dxx;
  const double dGxydV1y = a1*Dyy;
  const double dGxydV2x = a2*Dxy + b2*Dxx;
  const double dGxydV2y = a2*Dyy;
  
  const double dGyxdV1x = a1*Dxy + b1*Dxx;
  const double dGyxdV1y = a1*Dyy;
  const double dGyxdV2x = a2*Dxy + b2*Dxx;
  const double dGyxdV2y = a2*Dyy;
  
  const double dGyydV1x = 2*b1*Dxy;
  const double dGyydV1y = 2*b1*Dyy;
  const double dGyydV2x = 2*b2*Dxy;
  const double dGyydV2y = 2*b2*Dyy;
  
  // Calculate forces per area in rotated system: chain rule as in appendix C of Krüger (chain rule applied in eq. (C.13), but for the energy density).
  // Only two nodes are needed, third one is calcualted from momentum conservation
  // Note: If you calculate the derivatives in a straightforward manner, you get 8 terms (done here). Krueger exploits the symmetry of the G-matrix,
  //		which results in 6 elements, but with an additional factor 2 for the xy-elements (see also above at the definition of dI2dGxy).
  double f1_rot[2] = {0., 0.};
  double f2_rot[2] = {0., 0.};
  f1_rot[0] = - (dEdI1 * dI1dGxx * dGxxdV1x) - (dEdI1 * dI1dGxy * dGxydV1x) - (dEdI1 * dI1dGyx * dGyxdV1x) - (dEdI1 * dI1dGyy * dGyydV1x)
  - (dEdI2 * dI2dGxx * dGxxdV1x) - (dEdI2 * dI2dGxy * dGxydV1x) - (dEdI2 * dI2dGyx * dGyxdV1x) - (dEdI2 * dI2dGyy * dGyydV1x);
  f1_rot[1] = - (dEdI1 * dI1dGxx * dGxxdV1y) - (dEdI1 * dI1dGxy * dGxydV1y) - (dEdI1 * dI1dGyx * dGyxdV1y) - (dEdI1 * dI1dGyy * dGyydV1y)
  - (dEdI2 * dI2dGxx * dGxxdV1y) - (dEdI2 * dI2dGxy * dGxydV1y) - (dEdI2 * dI2dGyx * dGyxdV1y) - (dEdI2 * dI2dGyy * dGyydV1y);
  f2_rot[0] = - (dEdI1 * dI1dGxx * dGxxdV2x) - (dEdI1 * dI1dGxy * dGxydV2x) - (dEdI1 * dI1dGyx * dGyxdV2x) - (dEdI1 * dI1dGyy * dGyydV2x)
  - (dEdI2 * dI2dGxx * dGxxdV2x) - (dEdI2 * dI2dGxy * dGxydV2x) - (dEdI2 * dI2dGyx * dGyxdV2x) - (dEdI2 * dI2dGyy * dGyydV2x);
  f2_rot[1] = - (dEdI1 * dI1dGxx * dGxxdV2y) - (dEdI1 * dI1dGxy * dGxydV2y) - (dEdI1 * dI1dGyx * dGyxdV2y) - (dEdI1 * dI1dGyy * dGyydV2y)
  - (dEdI2 * dI2dGxx * dGxxdV2y) - (dEdI2 * dI2dGxy * dGxydV2y) - (dEdI2 * dI2dGyx * dGyxdV2y) - (dEdI2 * dI2dGyy * dGyydV2y);
  
  // Multiply by undeformed area
  f1_rot[0] *= A0;
  f1_rot[1] *= A0;
  f2_rot[0] *= A0;
  f2_rot[1] *= A0;
 

  // ****************** Wolfgang's version ***********
 /*
  // Left here for checking, but should be identical to the version above
  const double i11 = 1.0;
  const double i12 = 1.0;
  const double i21 = Gyy;
  const double i22 = -Gyx;
  const double i23 = i22;
  const double i24 = Gxx;

  //For sake of better readability shorten the call for the triangle's constants:
  A0 = iaparams->p.stretching_force_ibm.Area0;
  a1 = iaparams->p.stretching_force_ibm.a1;
  a2 = iaparams->p.stretching_force_ibm.a2;
  b1 = iaparams->p.stretching_force_ibm.b1;
  b2 = iaparams->p.stretching_force_ibm.b2;
  
  f1_rot[0] = A0*((-1)*e1*((i11*2*a1*dxx)+(i12*2*b1*dxy))+ (-1)*e2*((i21*2*a1*dxx)+(i22*(a1*dxy+b1*dxx))+(i23*(a1*dxy+b1*dxx))+(i24*2*b1*dxy)));
  f1_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b1*dyy))+ (-1)*e2*((i21*0.0)+(i22*a1*dyy)+(i23*a1*dyy)+(i24*2*b1*dyy)));
  
  f2_rot[0] = A0*((-1)*e1*((i11*2*a2*dxx)+(i12*2*b2*dxy))+ (-1)*e2*((i21*2*a2*dxx)+(i22*(a2*dxy+b2*dxx))+(i23*(a2*dxy+b2*dxx))+(i24*2*b2*dxy)));
  f2_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b2*dyy))+ (-1)*e2*((i21*0.0)+(i22*a2*dyy)+(i23*a2*dyy)+(i24*2*b2*dyy)));
  */
  
  // Rotate forces back into original position of triangle
  double force_1[3] = {0,0,0};
  double force_2[3] = {0,0,0};
  RotateForces(f1_rot,f2_rot, force_1,force_2, vec1, vec2);

  // set forces
  // forces are written in interface function via write_forces
  force = force_1;
  force2 = force_2;
  return 0;

}

//energy
int Bond::IbmTriel::calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
						     double *_energy) const {
  *_energy = 0.0;
  return 0;
}

// internal function
/** Rotate calculated trielastic forces in the 2d plane back to the 3d plane
 *Use knowledge that the x-axis in rotated system is parallel to r(p1->p2) in original system;
 *To find the corresponding unit vector to y in the rotated system, construct vector perpendicular to r(p1->p2);
 * note that f3 is not calculated here but is implicitly calculated by f3 = -(f1+f2) which is consistent with the literature
 */
void Bond::IbmTriel::RotateForces(const double f1_rot[2], const double f2_rot[2], double f1[3], double f2[3], double v12[3], double v13[3]) const {

  // fRot is in the rotated system, i.e. in a system where the side lPrime of the triangle (i.e. v12) is parallel to the x-axis, and the y-axis is perpendicular to the x-axis (cf. Krueger, Fig. 7.1c).
  //		I.e. fRot[XX] is the component parallel to the x-axis, fRot[YY]] the component parallel to the y-axis.
  // Now, the idea is to get these unit vectors for the x- and y-axis in the real coordinate system. They are named xu and yu below.
  // 		The x-component of the force in the real coordinate system is therefore: fRot[XX]*xu
  
  // xu is simple: The x-axis in the rotated system is parallel to v12 --> xu = v12 (+ normalization)
  double xu[3] = {0., 0., 0.};
  unit_vector(v12,xu);
  
  // yu needs to be orthogonal to xu, and point in the direction of node 3 in Krueger, Fig. 7.1b. Therefore:
  //		First get the projection of v13 onto v12: The direction is definied by xu, the length by the scalar product (scalar product can be interpreted as a projection, after all).
  //					--> sca * xu
  //		Then: v13 - sca * xu gives the component of v13 orthogonal to v12, i..e. perpendicular to the x-axis --> yu
  //		Last: Normalize yu.

  const double sca = scalar(v13,xu);
  double y[3] = {0., 0., 0.};
  for ( int i = 0; i < 3; i++)
    y[i] = v13[i] - (sca * xu[i]);

  double yu[3] = {0., 0., 0.};
  unit_vector(y,yu);
  
  // Calculate forces in 3D
  f1[0] = (f1_rot[0] * xu[0]) + (f1_rot[1] * yu[0]);
  f1[1] = (f1_rot[0] * xu[1]) + (f1_rot[1] * yu[1]);
  f1[2] = (f1_rot[0] * xu[2]) + (f1_rot[1] * yu[2]);
  
  f2[0] = (f2_rot[0] * xu[0]) + (f2_rot[1] * yu[0]);
  f2[1] = (f2_rot[0] * xu[1]) + (f2_rot[1] * yu[1]);
  f2[2] = (f2_rot[0] * xu[2]) + (f2_rot[1] * yu[2]);

}

int Bond::IbmTriel::ResetParams(double k1, double l0){
 
// Check if k1 is correct
  if ( fabs( m_k1 - k1) > 1e-9 ) {
    printf("k1 does not match while reading triel checkpoint!\n");
    return ES_ERROR;
  }
  
  // Check if l0 is correct
  if ( fabs( m_l0 - l0) > 1e-9 ) {
    printf("l0 does not match while reading triel checkpoint!\n");
    return ES_ERROR; }

  // Compute cache values a1, a2, b1, b2
  const double area0 = m_area0;
  const double lp0 = m_lp0;
  const double sinPhi0 = m_sinPhi0;
  const double cosPhi0 = m_cosPhi0;
  const double area2 = 2.0 * area0;
  const double a1 = -(l0*sinPhi0)/area2;
  const double a2 = - a1;
  const double b1 = (l0*cosPhi0 - lp0)/area2;
  const double b2 = -(l0*cosPhi0)/area2;
 
  // overide parameters
  m_a1 = a1;
  m_a2 = a2;
  m_b1 = b1;
  m_b2 = b2;

  return 0;

}
