/*
  Copyright (C) 2010 The ESPResSo project
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
#ifndef soft_H
#define soft_H

/** \file soft_sphere.h
 *  Routines to calculate the soft-sphere energy and/or  force 
 *  for a particle pair.
 *  \ref forces.c
*/

#ifdef SOFT_SPHERE

MDINLINE int printsoftIAToResult(Tcl_Interp *interp, int i, int j)
{
  char buffer[TCL_DOUBLE_SPACE];
  IA_parameters *data = get_ia_param(i, j);

  Tcl_PrintDouble(interp, data->soft_a, buffer);
  Tcl_AppendResult(interp, "soft-sphere ", buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_n, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_cut, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  Tcl_PrintDouble(interp, data->soft_offset, buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
  return TCL_OK;
}



/**Resultant Force due to a soft-sphere potential between two particles at interatomic separation r */
MDINLINE double soft_force_r(double a, double n, double r )
{
   return (a*n/pow(r, n+1));
}
/**Potential Energy due to a soft-sphere potential between two particles at interatomic separation r */
MDINLINE double soft_energy_r(double a, double n, double r )
{
   return (a/pow(r, n));
}


MDINLINE int soft_sphere_set_params(int part_type_a, int part_type_b,
                                    double a, double n, double cut, double offset)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* soft-sphere parameter should be symmetrical */
  data->soft_a      = data_sym->soft_a      = a;
  data->soft_n      = data_sym->soft_n      = n;
  data->soft_cut    = data_sym->soft_cut    = cut;
  data->soft_offset = data_sym->soft_offset = offset;
 
  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);


  return TCL_OK;
}



MDINLINE int soft_parser(Tcl_Interp * interp,
		       int part_type_a, int part_type_b,
		       int argc, char ** argv)
{
  /* parameters needed for soft-shere */
  double a, n, cut, offset;
  int change;

  /* get soft-sphere interaction type */
  if (argc < 5) {
    Tcl_AppendResult(interp, "soft-sphere potential needs 4 parameters: "
		     "<soft_a> <soft_n> <soft_cut> <soft_offset>",
		     (char *) NULL);
    return 0;
  }

  /* copy soft-sphere parameters */
  if ((! ARG_IS_D(1, a))     ||
      (! ARG_IS_D(2, n))     ||
      (! ARG_IS_D(3, cut))   ||
      (! ARG_IS_D(4, offset)   )) {
    Tcl_AppendResult(interp, "soft-sphere potential needs 4 parameters: "
		     "<soft_a> <soft_n> <soft_cut> <soft_offset>",
		     (char *) NULL);
    return TCL_ERROR;
  }
  change = 5;
	
  
  Tcl_ResetResult(interp);
  if (soft_sphere_set_params(part_type_a, part_type_b,
                             a, n, cut, offset) == TCL_ERROR) {
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  }
  return change;
}


/** Calculate soft-sphere potential force between particle p1 and p2 */
MDINLINE void add_soft_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist, double force[3])
{
  int j;
  double r_off, fac=0.0;
  if(dist < ia_params->soft_cut+ia_params->soft_offset) { 
    /* normal case: resulting force/energy smaller than zero. */
    r_off = dist - ia_params->soft_offset;
    if(r_off > 0.0) {
      fac = soft_force_r(ia_params->soft_a, ia_params->soft_n, r_off)/dist;
      for(j=0;j<3;j++)
	force[j] += fac * d[j];

#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: Soft_Sphere-Warning: Pair (%d-%d) force=%f dist=%f\n",
				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
    }
    

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: soft   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: soft   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));
  }
}

/** calculate soft-sphere energy between particle p1 and p2. */
MDINLINE double soft_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{
  double r_off;

  if(dist < ia_params->soft_cut+ia_params->soft_offset) {
    r_off = dist - ia_params->soft_offset;
    /* normal case: resulting force/energy smaller than zero. */
   
    return soft_energy_r(ia_params->soft_a, ia_params->soft_n, r_off);
    
  }
  return 0.0;
}



#endif /* ifdef SOFT_SPHERE */
#endif
