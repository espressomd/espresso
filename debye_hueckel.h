// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while u tilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef DEBYE_HUECKEL_H
#define DEBYE_HUECKEL_H
/** \file debye_hueckel.h
 *  Routines to calculate the Debye_Hueckel  Energy or/and Debye_Hueckel force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/

#ifdef ELECTROSTATICS

/** Structure to hold Debye-Hueckel Parameters. */
typedef struct {
  /** Cutoff for Debey-Hueckel interaction. */
  double r_cut;
  /** Debye kappa (inverse Debye length) . */
  double kappa;
} Debye_hueckel_params;

/** Structure containing the Debye-Hueckel parameters. */
extern Debye_hueckel_params dh_params;

/** \name Functions */
/************************************************************/
/*@{*/

MDINLINE int printdhToResult(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, dh_params.kappa, buffer);
  Tcl_AppendResult(interp, "dh ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, dh_params.r_cut, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

MDINLINE int dh_set_params(double kappa, double r_cut)
{
  if(dh_params.kappa < 0.0)
    return -1;

  if(dh_params.r_cut < 0.0)
    return -2;

  dh_params.kappa = kappa;
  dh_params.r_cut = r_cut;

  mpi_bcast_coulomb_params();

  return 1;
}

MDINLINE int inter_parse_dh(Tcl_Interp * interp, int argc, char ** argv)
{
  double kappa, r_cut;
  int i;

  if(argc < 2) {
    Tcl_AppendResult(interp, "Not enough parameters: inter coulomb dh <kappa> <r_cut>", (char *) NULL);
    return TCL_ERROR;
  }
  
  coulomb.method = COULOMB_DH;

  if(! ARG0_IS_D(kappa))
    return TCL_ERROR;
  if(! ARG1_IS_D(r_cut))
    return TCL_ERROR;

  if ( (i = dh_set_params(kappa, r_cut)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "dh kappa must be positiv.",(char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "dh r_cut must be positiv.",(char *) NULL);
      break;
    default:
      Tcl_AppendResult(interp, "unspecified error",(char *) NULL);
    }
    
    return TCL_ERROR;
  }

  return TCL_OK;
}

/** Computes the Debye_Hueckel pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
MDINLINE void add_dh_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
  int j;
  double kappa_dist, fac;
  
  if(dist < dh_params.r_cut) {
    if(dh_params.kappa > 0.0) {
      /* debye hueckel case: */
      kappa_dist = dh_params.kappa*dist;
      fac = coulomb.prefactor * p1->p.q * p2->p.q * (exp(-kappa_dist)/(dist*dist*dist)) * (1.0 + kappa_dist);
    }
    else {
      /* pure coulomb case: */
      fac = coulomb.prefactor * p1->p.q * p2->p.q / (dist*dist*dist);
    }
    for(j=0;j<3;j++)
      force[j] += fac * d[j];

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: DH   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: DH   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

MDINLINE double dh_coulomb_pair_energy(Particle *p1, Particle *p2, double dist)
{
  if(dist < dh_params.r_cut) {
    if(dh_params.kappa > 0.0)
      return coulomb.prefactor * p1->p.q * p2->p.q * exp(-dh_params.kappa*dist) / dist;
    else 
      return coulomb.prefactor * p1->p.q * p2->p.q / dist;
  }
  return 0.0;
}

/*@}*/
#endif

#endif
