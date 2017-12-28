#include "TabulatedBondDihedral.hpp"
#include "grid.hpp" //get_mi_vector
#include "utils.hpp" //vector_product

/** Calculate a tabulated dihedral force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1. p2, p3 and p4 and
    add it to the particle forces. This function is not tested yet.
    Needs feature TABULATED compiled in (see \ref config.hpp). */
    //force calculation
int Bond::TabulatedBondDihedral::calc_bonded_four_particle_force(Particle *p1, Particle *p2, 
								Particle *p3, Particle *p4, 
								double force1[3], double force2[3], 
								double force3[3], double force4[3]) 
  const

{
  int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi;
  /* force factors */
  double fac, f1[3], f4[3];

  /* dihedral angle */
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if (phi == -1.0) {
    for (i = 0; i < 3; i++) {
      force1[i] = 0.0;
      force2[i] = 0.0;
      force3[i] = 0.0;
    }
    return 0;
  }

  /* calculate force components (directions) */
  for (i = 0; i < 3; i++) {
    f1[i] = (v23Xv34[i] - cosphi * v12Xv23[i]) / l_v12Xv23;
    ;
    f4[i] = (v12Xv23[i] - cosphi * v23Xv34[i]) / l_v23Xv34;
  }
  vector_product(v23, f1, v23Xf1);
  vector_product(v23, f4, v23Xf4);
  vector_product(v34, f4, v34Xf4);
  vector_product(v12, f1, v12Xf1);

  /* table lookup */
  fac = m_tab_pot.force(phi);

  /* store dihedral forces */
  for (i = 0; i < 3; i++) {
    force1[i] = fac * v23Xf1[i];
    force2[i] = fac * (v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
    force3[i] = fac * (v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }

  return 0;
}

/** Calculate and return a tabulated dihedral energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1. p2,
    p3 and p4. This function is not tested yet. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
//energy calculation
int Bond::TabulatedBondDihedral::calc_bonded_four_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
				    Particle *p4, double *_energy) const
{
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);

  *_energy = m_tab_pot.energy(phi);

  return 0;
}

/** Calculates the dihedral angle between particle quadruple p1, p2,
p3 and p4. The dihedral angle is the angle between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). 
Vectors a, b and c are the bond vectors between consequtive particles.
If the a,b or b,c are parallel the dihedral angle is not defined in which
case the routine returns phi=-1. Calling functions should check for that
(Written by: Arijit Maitra) */
void Bond::TabulatedBondDihedral::calc_dihedral_angle(Particle *p1, Particle *p2, 
						      Particle *p3, Particle *p4, 
						      double a[3], double b[3], double c[3], 
						      double aXb[3], double *l_aXb, double bXc[3], 
						      double *l_bXc, double *cosphi, double *phi) 
  const
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
