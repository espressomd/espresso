// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
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


/** Calculates the dihedral angle between particel quadriple p1, p2,
p3 and p4. The dihedral angle is the angel between the planes
specified by the particle triples (p1,p2,p3) and (p2,p3,p4). Additional information is stored in:
 plane_vec1, plane_vec2, cos_v12_v23, cos_v23_v34*/
MDINLINE double calc_dihedral_angle(Particle *p1, Particle *p2, Particle *p3, Particle *p4, 
				    double plane_vec1[3], double *length_plane_vec1,
				    double plane_vec2[3], double *length_plane_vec2,
				    double *cos_v12_v23, double *cos_v23_v34, double *cosphi)
{
  int i;
  double tmp, tmpvec[3];
  /* vectors between particles */
  double v12[3], v23[3], v34[3];
  /* dihedral angle */
  double phi;

  get_mi_vector(v12, p2->r.p, p1->r.p);
  get_mi_vector(v23, p3->r.p, p2->r.p);
  get_mi_vector(v34, p4->r.p, p3->r.p);
  //  fprintf(stderr,"calc_dihedral_cos: v23=(%f,%f,%f)\n",v23[0],v23[1],v23[2]);
  //  fprintf(stderr,"calc_dihedral_cos: v34=(%f,%f,%f)\n",v34[0],v34[1],v34[2]);

  /* calculate plane_vec1 and plane_vec2 */
  tmp = 1.0/sqrlen(v23);
  *cos_v12_v23 = scalar(v12,v23)*tmp;
  *cos_v23_v34 = scalar(v23,v34)*tmp;
  //  fprintf(stderr,"calc_dihedral_cos: Angels (v12,v23)=%f, (v23,v34)=%f\n",
  //  acos(*cos_v12_v23),acos(*cos_v23_v34)); 
  for(i=0;i<3;i++) {
    plane_vec1[i] =  v12[i] - (*cos_v12_v23) * v23[i];
    plane_vec2[i] = -v34[i] + (*cos_v23_v34) * v23[i];
  }
  //fprintf(stderr,"plane_vec1=(%f,%f,%f)\n",plane_vec1[0],plane_vec1[1],plane_vec1[2]);
  //fprintf(stderr,"plane_vec2=(%f,%f,%f)\n",plane_vec2[0],plane_vec2[1],plane_vec2[2]);

  /* Calculate dihedral angle */
  *length_plane_vec1 = sqrt(sqrlen(plane_vec1));
  *length_plane_vec2 = sqrt(sqrlen(plane_vec2));
  if( *length_plane_vec1 == 0.0 && *length_plane_vec2 == 0.0 ) { return ANGLE_NOT_DEFINED; }
  /* return cosine of the dihedral angle */
  *cosphi = scalar(plane_vec1,plane_vec2)/((*length_plane_vec1)*(*length_plane_vec2));
  phi     = acos(*cosphi);

  /* take care of the degeneracy of the acos operation */
  vector_product(v12,v23,tmpvec);  
  if( scalar(tmpvec,v34) < 0.0 ) phi = (2.0*PI) - phi;

  return phi;
}


MDINLINE int calc_dihedral_force(Particle *p2, Particle *p1, Particle *p3, Particle *p4,
				 Bonded_ia_parameters *iaparams, double force2[3],
				 double force1[2], double force3[2])
{
  int i;
  /* vectors in plane of particle triples (p1,p2,p3) and (p2,p3,p4) and thier length. */
  double vec1[3], vecl1, vec2[3], vecl2;
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi, cos1223, cos2334;
  /* force factors */
  double fac, f1, f2, f4;

  phi = calc_dihedral_angle(p1, p2, p3, p4, vec1, &vecl1, vec2, &vecl2, 
			    &cos1223, &cos2334, &cosphi);
  /* fprintf(stderr,"add_dihedral_force: phi = %f, cosphi =%f\n",phi,cosphi); */
  
  if( phi == ANGLE_NOT_DEFINED ) {
    fprintf(stderr,"%d: calc_dihedral_force: undefined angle\n", this_node);
    return 1;
  }
  else {
    fac = - iaparams->p.dihedral.bend*iaparams->p.dihedral.phase;
    /* treat the different multiplicities */
    switch (iaparams->p.dihedral.mult) {
    case 1:
      break;
    case 2:
      fac *= 4.0*cosphi; break;
    case 3:
      fac *= 12.0*SQR(cosphi) - 3.0;  break;
    case 4:
      fac *= 16.0*cosphi * (2.0*SQR(cosphi) - 1.0 );  break;
    case 5:
      fac *= SQR(cosphi) * (80.0*SQR(cosphi) - 60.0) + 5.0;  break;
    case 6:
      fac *= ((SQR(cosphi) - 1.0) * 192.0*SQR(cosphi) + 36.0) * cosphi;  break; 
    }
    /* apply dihedral forces */
    vecl1 = 1.0/vecl1; vecl2 = 1.0/vecl2;
    for(i=0;i<3;i++) {
      f1 = fac * (vec2[i] - vec1[i]*cosphi) * vecl1;
      f4 = fac * (vec1[i] - vec2[i]*cosphi) * vecl2;
      f2 = (cos1223 - 1.0)*f1 - cos2334*f4;

      force1[i] = f1;
      force2[i] = f2;
      force3[i] = -(f1 + f2 + f4);
    }
  }
  return 0;
}

MDINLINE int dihedral_energy(Particle *p1, Particle *p2, Particle *p3, Particle *p4,
			     Bonded_ia_parameters *iaparams, double *_energy) 
{
  /* vectors in plane of particle triples (p1,p2,p3) and (p2,p3,p4) and thier length. */
  double vec1[3], vecl1, vec2[3], vecl2;
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi, cos1223, cos2334;
  /* force factors */
  double fac;

  phi = calc_dihedral_angle(p1, p2, p3, p4, vec1, &vecl1, vec2, &vecl2, 
			    &cos1223, &cos2334, &cosphi);

  if( phi != ANGLE_NOT_DEFINED ) {
    fac = iaparams->p.dihedral.phase;
    /* treat the different multiplicities */
    switch (iaparams->p.dihedral.mult) {
    case 1:
      fac *= cosphi; break;
    case 2:
      fac *= 2.0*SQR(cosphi) - 1.0; break;
    case 3:
      fac *= (4.0*SQR(cosphi) - 3.0) * cosphi;  break;
    case 4:
      fac *= (8.0*SQR(cosphi) - 8.0 ) * SQR(cosphi) + 1;  break;
    case 5:
      fac *= ((16.0*SQR(cosphi) - 20.0) * SQR(cosphi) + 5.0) * cosphi;  break;
    case 6:
      fac *= ((32.0*SQR(cosphi) - 48.0) * SQR(cosphi) + 18.0) * SQR(cosphi) - 1;  break; 
    }
    fac += 1.0;
    *_energy = fac * iaparams->p.dihedral.bend;
    return 0;
  }
  /* angle undefined */
  fprintf(stderr,"%d: calc_dihedral_force: undefined angle\n", this_node);
  return 1;
}

#endif
