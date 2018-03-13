#include "Dihedral.hpp"
#include "grid.hpp" //for get_mi_vector
#include "utils.hpp" // for vector_product 
#include "interaction_data.hpp"

// for dihedral forces has to be added in a different way
void Bond::Dihedral::write_force_to_particle(Particle *p1, Particle *p2, Particle *p3, 
			     Particle *p4, double force[3], double force2[3],
			     double force3[3], double force4[3]) const 
{
  for (int j = 0; j < 3; j++) {
    p1->f.f[j] += force[j];
    p2->f.f[j] += force2[j];
    p3->f.f[j] += force3[j];
    p4->f.f[j] -= force[j] + force2[j] + force3[j];
  };

}


/** Calculates the dihedral angle between particle quadruple p1, p2,
p3 and p4. The dihedral angle is the angle between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). 
Vectors a, b and c are the bond vectors between consequtive particles.
If the a,b or b,c are parallel the dihedral angle is not defined in which
case the routine returns phi=-1. Calling functions should check for that
(Written by: Arijit Maitra) */
void Bond::Dihedral::calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
				  double a[3], double b[3], double c[3], 
				  double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
				  double *cosphi, double *phi) const
{
  int i;

  get_mi_vector(a, p2->r.p, p1->r.p);
  get_mi_vector(b, p3->r.p, p2->r.p);
  get_mi_vector(c, p4->r.p, p3->r.p);

  /* calculate vector product a X b and b X c */
  vector_product(a, b, aXb);
  vector_product(b, c, bXc);

  /* calculate the unit vectors */
  *l_aXb = sqrt(sqrlen(aXb));
  *l_bXc = sqrt(sqrlen(bXc));

  /* catch case of undefined dihedral angle */
   if ( *l_aXb <= TINY_LENGTH_VALUE || *l_bXc <= TINY_LENGTH_VALUE ) { *phi = -1.0; *cosphi = 0; return;}

  for (i=0;i<3;i++) {
    aXb[i] /= *l_aXb;
    bXc[i] /= *l_bXc;
  }

  *cosphi = scalar(aXb, bXc);

  if ( fabs(fabs(*cosphi)-1)  < TINY_SIN_VALUE  ) *cosphi = dround(*cosphi);

  /* Calculate dihedral angle */
  *phi = acos(*cosphi);
  if( scalar(aXb, c) < 0.0 ) *phi = (2.0*PI) - *phi;

}


//force calculation
/** calculate dihedral force between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo,
    more general new dihedral form by Ana.
*/
//p1->p2, p2->p1, force->force2, force2-> force
int Bond::Dihedral::calc_bonded_four_particle_force(Particle *p2, Particle *p1, Particle *p3, 
				   Particle *p4, double force2[3], double force[3], 
						   double force3[3], double force4[3]) const {

  int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi, sinmphi_sinphi;
  /* force factors */
  double fac, f1[3], f4[3];

  /* dihedral angle */
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if ( phi == -1.0 ) { 
    for(i=0;i<3;i++) { force[i] = 0.0; force2[i] = 0.0; force3[i] = 0.0; }
    return 0;
  }

  /* calculate force components (directions) */
  for(i=0;i<3;i++)  {
    f1[i] = (v23Xv34[i] - cosphi*v12Xv23[i])/l_v12Xv23;;
    f4[i] = (v12Xv23[i] - cosphi*v23Xv34[i])/l_v23Xv34;
  }
  vector_product(v23, f1, v23Xf1);
  vector_product(v23, f4, v23Xf4);
  vector_product(v34, f4, v34Xf4);
  vector_product(v12, f1, v12Xf1);

  /* calculate force magnitude */
#ifdef OLD_DIHEDRAL
  fac = m_bend * m_phase * m_mult;
#else
  fac = -m_bend * m_mult;
#endif

  if(fabs(sin(phi)) < TINY_SIN_VALUE) {
#ifdef OLD_DIHEDRAL
    sinmphi_sinphi = m_mult * cos(2.0*PI - m_mult*phi)/cos(phi);
#else
    /*(comes from taking the first term of the MacLaurin expansion of
      sin(n*phi - phi0) and sin(phi) and then making the division).
      The original code had a 2PI term in the cosine (cos(2PI - nPhi))
      but I removed it because it wasn't doing anything. AnaVV*/
    sinmphi_sinphi = m_mult*
      cos(m_mult*phi - m_phase)/cosphi;
#endif
  }
  else {
#ifdef OLD_DIHEDRAL
    sinmphi_sinphi = sin(m_mult*phi)/sin(phi);
#else
    sinmphi_sinphi = sin(m_mult*phi - m_phase)/sin(phi);
#endif
  }

  fac *= sinmphi_sinphi;

  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force[i] = fac*v23Xf1[i];
      force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }
  return 0;
}

//energy calculation
//p1->p2, p2->p1
/** calculate dihedral energy between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo */
int Bond::Dihedral::calc_bonded_four_particle_energy(Particle *p2, Particle *p1, Particle *p3, 
						    Particle *p4, double *_energy) const {

  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;
  /* energy factors */
  double fac;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
#ifdef OLD_DIHEDRAL
  fac = m_phase * cos(m_mult*phi);
#else
  fac = -cos(m_mult*phi -m_phase);
#endif
  fac += 1.0;
  fac *=  m_bend;

  *_energy = fac;

  return 0;

}

boost::any Bond::Dihedral::get_bond_parameters_from_bond() const
{

  Dihedral_bond_parameters params = {m_mult, m_bend, m_phase};
  return boost::any(params);
  
}
