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
  /** ionic strength . */
  double kappa;
  /** epsilon1 (continuum dielectric constant inside) . */
  double epsilon1;
  /** epsilon2 (continuum dielectric constant outside) . */
  double epsilon2;
  /** Cutoff for Reaction Field interaction. */
  double r_cut;
  /** B important prefactor . */
  double B;
} Reaction_field_params;

/** Structure containing the Reaction Field parameters. */
extern Reaction_field_params rf_params;

/** \name Functions */
/************************************************************/
/*@{*/

MDINLINE int printrfToResult(Tcl_Interp *interp)
{
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, rf_params.kappa, buffer);
  Tcl_AppendResult(interp, "rf ", buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, rf_params.epsilon1, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, rf_params.epsilon2, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, rf_params.r_cut, buffer);
  Tcl_AppendResult(interp, buffer, (char *) NULL);

  return TCL_OK;
}

MDINLINE int rf_set_params(double kappa,double epsilon1,double epsilon2, double r_cut)
{
  rf_params.kappa = kappa;
  rf_params.epsilon1 = epsilon1;
  rf_params.epsilon2 = epsilon2;
  rf_params.r_cut = r_cut;
  rf_params.B =(2*(epsilon1-epsilon2)*(1+kappa*r_cut)-epsilon2*kappa*kappa*r_cut*r_cut)/((epsilon1+2*epsilon2)*(1+kappa*r_cut)+epsilon2*kappa*kappa*r_cut*r_cut);
  if(rf_params.epsilon1 < 0.0)
    return -1;

  if(rf_params.epsilon2 < 0.0)
    return -1;

  if(rf_params.r_cut < 0.0)
    return -2;

  mpi_bcast_coulomb_params();

  return 1;
}

MDINLINE int inter_parse_rf(Tcl_Interp * interp, int argc, char ** argv)
{
  double kappa,epsilon1,epsilon2, r_cut;
  int i;

  if(argc < 4) {
    Tcl_AppendResult(interp, "rf needs 4 parameters: "
                               "<kappa> <epsilon1> <epsilon2> <r_cut>",(char *) NULL);
    return TCL_ERROR;
  }

  coulomb.method = COULOMB_RF;

  if ((! ARG_IS_D(0, kappa))      ||
      (! ARG_IS_D(1, epsilon1))   ||
      (! ARG_IS_D(2, epsilon2))   ||
      (! ARG_IS_D(3, r_cut)        )) {
      Tcl_AppendResult(interp, "rf needs 4 parameters: "
                               "<kappa> <epsilon1> <epsilon2> <r_cut>",(char *) NULL);
       return TCL_ERROR;
  }

  if ( (i = rf_set_params(kappa,epsilon1,epsilon2,r_cut)) < 0) {
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
  double fac;

#ifdef WATER
  //Change
  double p1_com[3],p2_com[3],com_dist;

#ifndef WATER_FLEX
  if (p1->p.mol_id==p2->p.mol_id) return;
#endif

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }

  if (com_dist < rf_params.r_cut){
#else
  if (dist < rf_params.r_cut) {
#endif
     /*reaction field prefactor*/
     fac = 1.0 / (dist*dist*dist)  +  rf_params.B / (rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
     fac *= coulomb.prefactor * p1->p.q * p2->p.q;

     for (j=0;j<3;j++)
         force[j] += fac * d[j];

     ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
     ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

MDINLINE double rf_coulomb_pair_energy(Particle *p1, Particle *p2, double dist)
{
  double  fac;
#ifdef WATER
  //Change H2O
  double p1_com[3],p2_com[3],com_dist;

#ifndef WATER_FLEX
  if (p1->p.mol_id==p2->p.mol_id) return 0.0;
#endif

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return 0.0;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }
  if (com_dist < rf_params.r_cut) {
#else
  if (dist < rf_params.r_cut) {
#endif
     fac = 1.0 / dist  -  (rf_params.B*dist*dist) / (2*rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
     //cut off part
     fac -= (1-rf_params.B/2)  / rf_params.r_cut;
     fac *= coulomb.prefactor * p1->p.q * p2->p.q;
     return fac;
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
  double B;
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  B=(2*(epsilon1-epsilon2)*(1+kappa*r_cut)-epsilon2*kappa*kappa*r_cut*r_cut)/((epsilon1+2*epsilon2)*(1+kappa*r_cut)+epsilon2*kappa*kappa*r_cut*r_cut);

  /* LJ should be symmetrically */
  data->rf_coul_pref = data_sym->rf_coul_pref = coul_pref;
  data->rf_kappa     = data_sym->rf_kappa     = kappa;
  data->rf_epsilon1  = data_sym->rf_epsilon1  = epsilon1;
  data->rf_epsilon2  = data_sym->rf_epsilon2  = epsilon2;
  data->rf_r_cut     = data_sym->rf_r_cut     = r_cut;
  data->rf_B         = data_sym->rf_B         = B;

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
#ifdef WATER
  //Change
  double p1_com[3],p2_com[3],com_dist;

#ifndef WATER_FLEX
  if (p1->p.mol_id==p2->p.mol_id) return;
#endif

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }
  if (com_dist < ia_params->rf_r_cut) {
#else
  if (dist < ia_params->rf_r_cut) {
#endif
     /*reaction field prefactor*/
     fac = 1.0 / (dist*dist*dist)  +  ia_params->rf_B / (ia_params->rf_r_cut*ia_params->rf_r_cut*ia_params->rf_r_cut);
     fac *= ia_params->rf_coul_pref * p1->p.q * p2->p.q;
     for (j=0;j<3;j++)
         force[j] += fac * d[j];

     ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
     ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

MDINLINE double interrf_pair_energy(Particle *p1, Particle *p2,IA_parameters *ia_params, double dist)
{
  double fac;
#ifdef WATER
  //Change H2O
  double p1_com[3],p2_com[3],com_dist;

#ifndef WATER_FLEX
  if (p1->p.mol_id==p2->p.mol_id) return 0.0;
#endif

  if ((get_com_h2o(p1,p1_com) == -1 ) || (get_com_h2o(p2,p2_com)==-1)){
     return 0.0;
  }
  else{
     com_dist=min_distance(p1_com,p2_com);
  }

  if (com_dist < ia_params->rf_r_cut) {
#else
  if (dist < ia_params->rf_r_cut) {
#endif
     fac = 1.0 / dist  -  (ia_params->rf_B*dist*dist) / (2*ia_params->rf_r_cut*ia_params->rf_r_cut*ia_params->rf_r_cut);
     //cut off part
     fac -= (1-ia_params->rf_B/2)  / ia_params->rf_r_cut;
     fac *= ia_params->rf_coul_pref * p1->p.q * p2->p.q;
     return fac;
  }
  return 0.0;
}

#endif

/*@}*/
#endif

#endif
