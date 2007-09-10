// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while u tilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef REACTION_FIELD_H
#define REACTION_FIELD_H
/** \file reaction_field.h
 *  Routines to calculate the Reaction Field Energy or/and force 
 *  for a particle pair.
 *  M. Neumann, J. Chem. Phys 82, 5663 (1985)
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:junghans@mpip-mainz.mpg.de">Christoph</a>
    stolen form Matej's ACG version of Espresso
*/

#ifdef ELECTROSTATICS

#include "statistics.h"

/** Structure to hold Reaction Field Parameters. */
typedef struct {
  /** Cutoff for Reaction Field interaction. */
  double r_cut;
  /** eps (continuum dielectric constant) . */
  double eps;
} Reaction_field_params;

/** Structure containing the Reaction Field parameters. */
extern Reaction_field_params rf_params;

/** \name Functions */
/************************************************************/
/*@{*/

MDINLINE int printrfToResult(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, rf_params.eps, buffer);
  Tcl_AppendResult(interp, "rf ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, rf_params.r_cut, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

MDINLINE int rf_set_params(double eps, double r_cut)
{
  if(rf_params.eps < 0.0)
    return -1;

  if(rf_params.r_cut < 0.0)
    return -2;

  rf_params.eps = eps;
  rf_params.r_cut = r_cut;

  mpi_bcast_coulomb_params();

  return 1;
}

MDINLINE int inter_parse_rf(Tcl_Interp * interp, int argc, char ** argv)
{
  double eps, r_cut;
  int i;

  if(argc < 2) {
    Tcl_AppendResult(interp, "Not enough parameters: inter coulomb rf <eps> <r_cut>", (char *) NULL);
    return TCL_ERROR;
  }
  
  coulomb.method = COULOMB_RF;

  if(! ARG0_IS_D(eps))
    return TCL_ERROR;
  if(! ARG1_IS_D(r_cut))
    return TCL_ERROR;

  if ( (i = rf_set_params(eps, r_cut)) < 0) {
    switch (i) {
    case -1:
      Tcl_AppendResult(interp, "rf eps must be positive.",(char *) NULL);
      break;
    case -2:
      Tcl_AppendResult(interp, "rf r_cut must be positive.",(char *) NULL);
      break;
    default:
      Tcl_AppendResult(interp, "unspecified error",(char *) NULL);
    }
    
    return TCL_ERROR;
  }

  return TCL_OK;
}

/** Computes the Reaction Field pair force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param d         Vector pointing from p1 to p2.
    @param dist      Distance between p1 and p2.
    @param force     returns the force on particle 1.
*/
MDINLINE void add_rf_coulomb_pair_force(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
  int j;
  double B0, fac, fac_rf;

  if(dist < rf_params.r_cut) {
    /*reaction field prefactor*/
    B0 = 2.0*(rf_params.eps-1.0)/(2.0*rf_params.eps+1.0);
    fac_rf = -B0*coulomb.prefactor * p1->p.q * p2->p.q / (rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
    fac = coulomb.prefactor * p1->p.q * p2->p.q / (dist*dist*dist);

    for(j=0;j<3;j++)
       force[j] += (fac+fac_rf) * d[j];

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}






MDINLINE double rf_coulomb_pair_energy(Particle *p1, Particle *p2, double dist)
{
  double  B0, fac, fac_cut,fac_rf;
  if(dist < rf_params.r_cut) {
    B0 = 2.0*(rf_params.eps-1.0)/(2.0*rf_params.eps+1.0);
    //no minus here
    fac_rf = B0*coulomb.prefactor * p1->p.q * p2->p.q * dist * dist   / (2.0 * rf_params.r_cut * rf_params.r_cut * rf_params.r_cut); 
    // also no minus
    fac =       coulomb.prefactor * p1->p.q * p2->p.q              / dist;
    // cutoff part -- fac+fac_rf at dist=r_cut
    fac_cut = -  coulomb.prefactor *p1->p.q * p2->p.q *(1.0+B0/2.0)/rf_params.r_cut;
    return fac+fac_rf+fac_cut;
  }
  return 0.0;
}

/*@}*/
#endif

#endif
