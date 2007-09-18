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

#ifdef WATER
  //Change
  double p1_com[3],p2_com[3],com_dist;

  if (p1->p.mol_id==p2->p.mol_id) return;

  get_com_h2o(p1,p1_com);
  get_com_h2o(p2,p2_com);
  com_dist=min_distance(p1_com,p2_com);

  if(com_dist < rf_params.r_cut) {
#else
  if(dist < rf_params.r_cut) {
#endif
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

#ifdef WATER
  //Change H2O
  double p1_com[3],p2_com[3],com_dist;

  if (p1->p.mol_id==p2->p.mol_id) return 0.0;

  get_com_h2o(p1,p1_com);
  get_com_h2o(p2,p2_com);
  com_dist=min_distance(p1_com,p2_com);

  if(com_dist < rf_params.r_cut) {
#else
  if(dist < rf_params.r_cut) {
#endif
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

#ifdef INTER_RF
/*from I. G. Tironi et al., J. Chem. Phys. 102, 5451 (1995)*/
MDINLINE int printinterrfIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->rf_coul_pref, buffer);
  Tcl_AppendResult(interp, "inter_rf ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->rf_kappa, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->rf_epsilon1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->rf_epsilon2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->rf_r_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);

  return TCL_OK;
}

MDINLINE int interrf_set_params(int part_type_a, int part_type_b,double coul_pref,
				      double kappa, double epsilon1, double epsilon2,
				      double r_cut)
{
  double B0,B1;
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  B0=(epsilon1-2*epsilon2*(1+kappa*r_cut))/(epsilon2*(1+kappa*r_cut));
  B1=((epsilon1-4*epsilon2)*(1+kappa*r_cut)-2*epsilon2*kappa*kappa*r_cut*r_cut)/((epsilon1+2*epsilon2)*(1+kappa*r_cut)+epsilon2*kappa*kappa*r_cut*r_cut);

  /* LJ should be symmetrically */
  data->rf_coul_pref = data_sym->rf_coul_pref = coul_pref;
  data->rf_kappa     = data_sym->rf_kappa     = kappa;
  data->rf_epsilon1  = data_sym->rf_epsilon1  = epsilon1;
  data->rf_epsilon2  = data_sym->rf_epsilon2  = epsilon2;
  data->rf_r_cut     = data_sym->rf_r_cut     = r_cut;
  data->rf_B0        = data_sym->rf_B0        = B0;
  data->rf_B1        = data_sym->rf_B1        = B1;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int interrf_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for RF */
  double coul_pref,kappa,epsilon1,epsilon2,r_cut;
  int change;

  /* get lennard-jones interaction type */
  if (argc < 6) {
    Tcl_AppendResult(interp, "inter_rf needs 5 parameters: "
		     "<coul_pref> <kappa> <epsilon1> <epsilon2> <r_cut>",
		     (char *) NULL);
    return 0;
  }

  /* copy lennard-jones parameters */
  if ((! ARG_IS_D(1, coul_pref))      ||
      (! ARG_IS_D(2, kappa))   ||
      (! ARG_IS_D(3, epsilon1))   ||
      (! ARG_IS_D(4, epsilon2))   ||
      (! ARG_IS_D(5, r_cut)        )) {
    Tcl_AppendResult(interp, "inter_rf needs 5 parameters: "
		     "<coul_pref> <kappa> <epsilon1> <epsilon2> <r_cut>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 6;
	
  if (interrf_set_params(part_type_a, part_type_b,coul_pref,
			       kappa,epsilon1,epsilon2,r_cut
			       ) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}
MDINLINE void add_interrf_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double fac;
  if(dist < ia_params->rf_r_cut) {
    /*reaction field prefactor*/
    fac = (1.0/(dist*dist*dist)-(1+ia_params->rf_B1)/(ia_params->rf_r_cut*ia_params->rf_r_cut*ia_params->rf_r_cut));
    fac *= ia_params->rf_coul_pref * p1->p.q * p2->p.q;
    for(j=0;j<3;j++)
       force[j] += fac * d[j];

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

MDINLINE double interrf_pair_energy(Particle *p1, Particle *p2,IA_parameters *ia_params, double dist)
{
  double fac;
  if(dist < ia_params->rf_r_cut) {
    fac =ia_params->rf_coul_pref * p1->p.q * p2->p.q*(1.0/dist+(1+ia_params->rf_B0)/ia_params->rf_r_cut);
    return fac;
  }
  return 0.0;
}

#endif

/*@}*/
#endif

#endif
