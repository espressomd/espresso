/*
  Copyright (C) 2010,2011 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef REACTION_FIELD_H
#define REACTION_FIELD_H
/** \file reaction_field.h
 *  Routines to calculate the Reaction Field Energy or/and force 
 *  for a particle pair.
 *  M. Neumann, J. Chem. Phys 82, 5663 (1985)
 *  \ref forces.c
 *
*/

#ifdef ELECTROSTATICS

#include "mol_cut.h"

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

MDINLINE int tclprint_to_result_rf(Tcl_Interp *interp,char *name)
{
  char buffer[TCL_DOUBLE_SPACE];
  sprintf(buffer,"%s",name);
  Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
  Tcl_PrintDouble(interp, rf_params.kappa, buffer);
  Tcl_AppendResult(interp, buffer, " ",(char *) NULL);
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

MDINLINE int tclcommand_inter_coulomb_parse_rf(Tcl_Interp * interp, int argc, char ** argv,int method)
{
  double kappa,epsilon1,epsilon2, r_cut;
  int i;

  if(argc < 4) {
    Tcl_AppendResult(interp, "rf needs 4 parameters: "
                               "<kappa> <epsilon1> <epsilon2> <r_cut>",(char *) NULL);
    return TCL_ERROR;
  }

  coulomb.method = method;

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

MDINLINE void add_rf_coulomb_pair_force_no_cutoff(Particle *p1, Particle *p2, double d[3], double dist, double force[3])
{
   int j;
   double fac;
   fac = 1.0 / (dist*dist*dist)  +  rf_params.B / (rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
   fac *= coulomb.prefactor * p1->p.q * p2->p.q;
   
   for (j=0;j<3;j++)
         force[j] += fac * d[j];
   
   ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
   ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
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
   if (dist < rf_params.r_cut)
   {
      add_rf_coulomb_pair_force_no_cutoff(p1,p2,d,dist,force);
   }
}

MDINLINE double rf_coulomb_pair_energy_no_cutoff(Particle *p1, Particle *p2, double dist)
{
   double  fac;
   fac = 1.0 / dist  -  (rf_params.B*dist*dist) / (2*rf_params.r_cut*rf_params.r_cut*rf_params.r_cut);
   //cut off part
   fac -= (1-rf_params.B/2)  / rf_params.r_cut;
   fac *= coulomb.prefactor * p1->p.q * p2->p.q;
   return fac;
}

MDINLINE double rf_coulomb_pair_energy(Particle *p1, Particle *p2, double dist)
{
  if (dist < rf_params.r_cut)
  {
     return rf_coulomb_pair_energy_no_cutoff(p1,p2,dist);
  }
  return 0.0;
}

/*from I. G. Tironi et al., J. Chem. Phys. 102, 5451 (1995)*/
#ifdef INTER_RF
MDINLINE int tclprint_to_result_interrfIA(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  sprintf(buffer,"%i",data->rf_on);
  Tcl_AppendResult(interp, "inter_rf ", buffer, " ",(char *) NULL);
  return TCL_OK;
}

MDINLINE int interrf_set_params(int part_type_a, int part_type_b,int rf_on)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  data->rf_on         = data_sym->rf_on         = rf_on;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  return TCL_OK;
}

MDINLINE int tclcommand_inter_parse_interrf(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for RF */
  int rf_on;
  int change;

  /* get reaction_field interaction type */
  if (argc < 2) {
    Tcl_AppendResult(interp, "inter_rf needs 1 parameter: "
		     "<rf_on>",
		     (char *) NULL);
    return 0;
  }

  /* copy reaction_field parameters */
  if (! ARG_IS_I(1, rf_on)) {
    Tcl_AppendResult(interp, "<rf_on> must be int",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 2;
	
  if (! ((rf_on==0) || (rf_on==1)) ) {
    Tcl_AppendResult(interp, "rf_on must be 0 or 1", (char *) NULL);
    return 0;
  }
  if (interrf_set_params(part_type_a, part_type_b,rf_on) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}

MDINLINE void add_interrf_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
#ifdef ONEPART_DEBUG
  double fac=0.0 ; /* TODO: this  variable was not declared. Now the code compiles, but I have no idea of what value to assign to it (MS) */
#endif
  if ((ia_params->rf_on ==1) && CUTOFF_CHECK(dist < rf_params.r_cut)) {
     add_rf_coulomb_pair_force_no_cutoff(p1,p2,d, dist,force);
  }

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: INTER_RF   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
}

MDINLINE double interrf_pair_energy(Particle *p1, Particle *p2,IA_parameters *ia_params, double dist)
{
  double val;
  if ((ia_params->rf_on==1) && CUTOFF_CHECK(dist < rf_params.r_cut)) {
     val=rf_coulomb_pair_energy_no_cutoff(p1,p2,dist);
     return val;
  }
  return 0.0;
}

#endif

/*@}*/
#endif

#endif
