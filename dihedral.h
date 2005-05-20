// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
#ifndef DIHEDRAL_H
#define DIHEDRAL_H
/** \file dihedral.h Routines to calculate the dihedral energy or/and
 *  and force for a particle quadruple.  Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/

#define ANGLE_NOT_DEFINED -100




/// set dihedral parameters
MDINLINE int dihedral_set_params(int bond_type, int mult, double bend, double phase)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_DIHEDRAL;
  bonded_ia_params[bond_type].num  = 3;
  bonded_ia_params[bond_type].p.dihedral.mult = mult;
  bonded_ia_params[bond_type].p.dihedral.bend = bend;
  bonded_ia_params[bond_type].p.dihedral.phase = phase;

  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/** Calculates the dihedral angle between particle quadruple p1, p2,
p3 and p4. The dihedral angle is the angle between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). 
Vectors a, b and c are the bond vectors between consequtive particles.
(Written by: Arijit Maitra) */
MDINLINE void calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
				  double a[3], double b[3], double c[3], 
				  double aXb[3], double *l_aXb, double bXc[3], double *l_bXc, 
				  double *cosphi, double *phi)
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

/** calculate dihedral force between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo */
MDINLINE int calc_dihedral_force(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force2[3],
				 double force1[2], double force3[2])
{
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
  fac =   iaparams->p.dihedral.bend * iaparams->p.dihedral.phase * iaparams->p.dihedral.mult;
  
  if(fabs(sin(phi)) < TINY_SIN_VALUE) {
    sinmphi_sinphi =  iaparams->p.dihedral.mult * cos(2.0*PI -  iaparams->p.dihedral.mult*phi)/cos(phi); 
  }
  else {
    sinmphi_sinphi = sin( iaparams->p.dihedral.mult*phi)/sin(phi);
  }

  fac *= sinmphi_sinphi;


  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force1[i] = fac*v23Xf1[i];
      force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }
  return 0;
}

/** calculate dihedral energy between particles p1, p2 p3 and p4 
    Written by Arijit Maitra, adapted to new force interface by Hanjo */
MDINLINE int dihedral_energy(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
			     Bonded_ia_parameters *iaparams, double *_energy) 
{
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;
  /* energy factors */
  double fac;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
  
  fac =  iaparams->p.dihedral.phase * cos( iaparams->p.dihedral.mult*phi);
  fac += 1.0;
  fac *=  iaparams->p.dihedral.bend;

  return fac;
}

#endif
