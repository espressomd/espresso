// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef BUCKINGHAM_H
#define BUCKINGHAM_H
/** \file buckingham.h
 *  Routines to calculate the Buckingham energy and/or  force
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:arijitmaitra@uni-muenster.de">Arijit</a>
*/

#ifdef BUCKINGHAM

MDINLINE int printbuckIAToResult(Tcl_Interp *interp, int i, int j)
{
    char buffer[TCL_DOUBLE_SPACE];
    IA_parameters *data = get_ia_param(i, j);

    Tcl_PrintDouble(interp, data->BUCK_A, buffer);
    Tcl_AppendResult(interp, "buckingham ", buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_B, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_C, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_D, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_cut, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_discont, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_shift, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_capradius, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_F1, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
    Tcl_PrintDouble(interp, data->BUCK_F2, buffer);
    Tcl_AppendResult(interp, buffer, " ", (char *) NULL);
 
    return TCL_OK;
}

MDINLINE int buckforcecap_set_params(double buckforcecap)
{
  if (buck_force_cap != -1.0)
    mpi_buck_cap_forces(buck_force_cap);

  return TCL_OK;
}


/**Resultant Force due to a buckingham potential between two particles at interatomic separation r greater than or equal to discont*/
MDINLINE double buck_force_r(double A, double B, double C, double D, double r )
{
   return (A*B*exp(-B*r) - 6.0*C/pow(r, 7) - 4.0*D/pow(r, 5));
}
/**Potential Energy due to a buckingham potential between two particles at interatomic separation r greater than or equal to discont*/
MDINLINE double buck_energy_r(double A, double B, double C, double D, double shift, double r )
{
   return (A*exp(-B*r) - C/pow(r, 6) - D/pow(r, 4) + shift);
}



MDINLINE int buckingham_set_params(int part_type_a, int part_type_b,
               			   double A, double B, double C, double D, double cut,
		         	   double discont, double shift, double cap_radius,
			           double F1, double F2)
{
  IA_parameters *data, *data_sym;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);

  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym) {
    return TCL_ERROR;
  }

  /* BUCKINGHAM should be symmetrical */
  data->BUCK_A      = data_sym->BUCK_A      = A;
  data->BUCK_B      = data_sym->BUCK_B      = B;
  data->BUCK_C      = data_sym->BUCK_C      = C;
  data->BUCK_D      = data_sym->BUCK_D      = D;
  data->BUCK_cut    = data_sym->BUCK_cut    = cut;
  data->BUCK_discont    = data_sym->BUCK_discont    = discont;
  data->BUCK_shift  = data_sym->BUCK_shift  = shift;

  if (cap_radius > 0.0) {
    data->BUCK_capradius = cap_radius;
    data_sym->BUCK_capradius = cap_radius;
  }
  /* Replace the buckingham potential for interatomic dist. less
    than or equal to discontinuity by a straight line (F1+F2*r) */
  if (F1==0.0 && F2==0.0) {
     F1 = buck_energy_r(A, B, C, D, shift, discont) +
          discont*buck_force_r(A, B, C, D, discont);
     F2 = -buck_force_r(A, B, C, D, discont);
  }
  data->BUCK_F1 = data_sym->BUCK_F1=F1;
  data->BUCK_F2 = data_sym->BUCK_F2=F2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (buck_force_cap != -1.0)
     mpi_buck_cap_forces(buck_force_cap);

  return TCL_OK;
}




///parser for the forcecap
MDINLINE int inter_parse_buckforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];

  if (argc == 0) {
    if (buck_force_cap == -1.0)
      Tcl_AppendResult(interp, "buckforcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, buck_force_cap, buffer);
      Tcl_AppendResult(interp, "buckforcecap ", buffer, (char *) NULL);
    }
    return TCL_OK;
  }

  if (argc > 1) {
    Tcl_AppendResult(interp, "inter buckforcecap takes at most 1 parameter",
		     (char *) NULL);
    return TCL_ERROR;
  }

  if (ARG0_IS_S("individual"))
      buck_force_cap = -1.0;
  else if (! ARG0_IS_D(buck_force_cap) || buck_force_cap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(buckforcecap_set_params(buck_force_cap),
	      "If you can read this, you should change it. (Use the source Luke!)");
  return TCL_ERROR;
}





MDINLINE int buckingham_parser(Tcl_Interp * interp,
                         int part_type_a, int part_type_b,
                         int argc, char ** argv)
{
  /* parameters needed for buckingham */
  double buck_A,buck_B,buck_C,buck_D,buck_cut,buck_discont,buck_shift,buck_cap_radius,F1,F2;
  int change;


  /* get buckingham interaction type */
  if (argc < 8) {
     Tcl_AppendResult(interp, "buckingham needs 7 parameters: "
		      "<buck_A> <buck_B> <buck_C> <buck_D> <buck_cut> <buck_discontinuity> <buck_shift> ",
		      (char *) NULL);
     return TCL_ERROR;
  }
  /* copy buckingham parameters */
  if ((! ARG_IS_D(1, buck_A))   ||
      (! ARG_IS_D(2, buck_B))   ||
      (! ARG_IS_D(3, buck_C))   ||
      (! ARG_IS_D(4, buck_D))   ||
      (! ARG_IS_D(5, buck_cut)) ||
      (! ARG_IS_D(6, buck_discont)) ||
      (! ARG_IS_D(7, buck_shift))) {
    Tcl_AppendResult(interp, "buckingham needs 7 DOUBLE parameters: "
		     "<buck_A> <buck_B> <buck_C> <buck_D> <buck_cut> <buck_discontinuity> <buck_shift> ",
		     (char *) NULL);
      return TCL_ERROR;
  }
  change = 8;
  buck_cap_radius = -1.0;
  F1 = 0.0;
  F2 = 0.0;
  /* check whether there are additional doubles, cap radius, F1 and F2*/
  if (argc >= 9 && ARG_IS_D(8, buck_cap_radius))
  {
    change++;
    if(argc >= 10 && ARG_IS_D(9, F1))
    {
       change++;
       if(argc >= 11 && ARG_IS_D(10, F1))
       change++;
    }
  }
  else
    Tcl_ResetResult(interp);
  if(buck_discont>buck_cut)
  {
     Tcl_AppendResult(interp, "ERROR: <buck_cut> must be greater than <buck_discontinuity>",
		      (char *) NULL);
     return TCL_ERROR;
  }
  if (buckingham_set_params(part_type_a, part_type_b,
	        	    buck_A, buck_B, buck_C, buck_D,
			    buck_cut, buck_discont, buck_shift,
			    buck_cap_radius, F1, F2) == TCL_ERROR) {
     Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
     return 0;
  }
  return change;
}  



/** Calculate Buckingham force between particle p1 and p2 and add
    it to their force. */
MDINLINE void add_buck_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
                                  double d[3], double dist, double force[3])
{
  int j;
  double fac;
  if(dist < ia_params->BUCK_cut ) {
    if(ia_params->BUCK_capradius==0.0)
    {
     /* case: resulting force/energy greater than discontinuity and
              less than cutoff (true buckingham region) */
     if(dist > ia_params->BUCK_discont) {
      fac = buck_force_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, dist )/dist;
      for(j=0;j<3;j++) {
          force[j] += fac * d[j];
      }
#ifdef LJ_WARN_WHEN_CLOSE
      if(fac*dist > 1000) fprintf(stderr,"%d: BUCK-Warning: Pair (%d-%d) force=%f dist=%f\n",
      				  this_node,p1->p.identity,p2->p.identity,fac*dist,dist);
#endif
     }
     else
     {
      /* resulting force/energy in the linear region*/
         fac = -ia_params->BUCK_F2/dist;
         for(j=0;j<3;j++) {
	     force[j] += fac * d[j];
         }
     }
    }
    /* if Buckingham potential is capped. */
    else
    {
      if(dist>ia_params->BUCK_capradius) {
       fac = buck_force_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, dist )/dist;
       for(j=0;j<3;j++) {
	   force[j] += fac * d[j];
       }
     }
     else
     {
       fac = buck_force_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, ia_params->BUCK_capradius)/dist;
       for(j=0;j<3;j++) {
	   force[j] += fac * d[j];
       }
     }
   }

    ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: BUCK   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
    ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: BUCK   f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

    BUCK_TRACE(fprintf(stderr,"%d: BUCK: Pair (%d-%d) dist=%.3f: force+-: (%.3e,%.3e,%.3e)\n",
		     this_node,p1->p.identity,p2->p.identity,dist,fac*d[0],fac*d[1],fac*d[2]));
  }
}

/** calculate Buckingham energy between particle p1 and p2. */
MDINLINE double buck_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				double d[3], double dist)
{

  if(dist < ia_params->BUCK_cut) {
    if(ia_params->BUCK_capradius==0.0)
    {
     /* case: resulting force/energy greater than discont and
              less than cutoff (true buckingham region) */
     if(dist > ia_params->BUCK_discont)
        return buck_energy_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, ia_params->BUCK_shift, dist);
     else
        /* resulting force/energy in the linear region*/
        return (ia_params->BUCK_F1+ia_params->BUCK_F2*dist);
    }
    else
    {
        if(dist>ia_params->BUCK_capradius)
	  /*true buckigham region*/
	  return buck_energy_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, ia_params->BUCK_shift, dist);
	else
	  /*capped region*/
          return buck_energy_r(ia_params->BUCK_A, ia_params->BUCK_B, ia_params->BUCK_C, ia_params->BUCK_D, ia_params->BUCK_shift, ia_params->BUCK_capradius);
    }
  }
  return 0.0;
}

/** calculate buck_capradius from buckingham force cap */
MDINLINE void calc_buck_cap_radii(double force_cap)
{
  int i,j,cnt=0;
  IA_parameters *params = NULL;
  double force=0.0, frac2, frac6, frac8;
  double r0,r1 = 0.0,diff,exp_term,C_R,D_R;
  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      if(force_cap>0.0 ) {
        force = -params->BUCK_F2;
        if (force_cap<force)
	{
	  /* Solving numerically using Newton Raphson Technique */
	  force = 0.0;
	  cnt = 0;
	  r1 = (params->BUCK_cut+params->BUCK_discont)/2.0; //First guess value
	  r0 = 0.0 ;
	  while(fabs(r1 - r0)>1.0e-10) {
	     r0 = r1;
	     exp_term = params->BUCK_A*params->BUCK_B*exp(-params->BUCK_B*r0);
	     frac2 = SQR(r0);
	     frac6 = frac2*frac2*frac2;
	     frac8 = frac6*frac2;
	     C_R = 6.0*params->BUCK_C/frac8;
	     D_R = 4.0*params->BUCK_D/frac6;
	     diff = (exp_term - C_R*r0 - D_R*r0 - force_cap)/(-params->BUCK_B*exp_term + 7.0*C_R + 5.0*D_R);
	     r1 = r0 - diff;
	     if(r1>params->BUCK_discont)
	        r1=0.5*(params->BUCK_discont+r0);
	     cnt++;
	     if(cnt>500)
	     {
	       fprintf(stderr,"%d: ERROR@buckingham.h: Failed to converge while determining Buckingham cap radius!!",this_node);
	       fprintf(stderr,"%d: tolerance = %f",this_node, diff);
               exit (0);
	     }
          }
	  frac2 = SQR(r1);
	  frac6 = frac2*frac2*frac2;
	  force = params->BUCK_A*params->BUCK_B*exp(-params->BUCK_B*r1) - 6.0*params->BUCK_C/(r1*frac6) - 4.0*params->BUCK_D*r1/(frac6);
          params->BUCK_capradius = r1;
	}
	else
	   params->BUCK_capradius = params->BUCK_discont;

      }
      else
	params->BUCK_capradius = 0.0;

      FORCE_TRACE(fprintf(stderr,"%d: Ptypes %d-%d have cap_radius %f and cap_force %f (iterations: %d)\n",
			  this_node,i,j,r1,force,cnt));
    }
  }
  BUCK_TRACE(fprintf(stderr,"%d: BUCK: Buckingham force cap imposed %f, Calculated force %f and Cap radius %f after %d iterations\n",this_node,force_cap,force,params->BUCK_capradius,cnt);
  );
}
#endif
#endif
