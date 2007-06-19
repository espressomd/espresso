// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef TAB_H
#define TAB_H

/** \file tab.h
 *  Routines to calculate the  energy and/or  force 
 *  for a particle pair or bonds via interpolating from lookup tables.
 *  \ref forces.c
 *  Needs feature TABULATED compiled in (see \ref config.h).

 *  <b>Responsible:</b>
 *  <a href="mailto:cooke@mpip-mainz.mpg.de">Ira</a>
*/

#ifdef TABULATED

#include "dihedral.h"

/// set parameters for the force capping of tabulated potentials
MDINLINE int tabforcecap_set_params(double tabforcecap)
{
  if (tab_force_cap != -1.0)
    mpi_tab_cap_forces(tab_force_cap);

  return TCL_OK;
}

/** Non-Bonded tabulated potentials:
    Reads tabulated parameters and force and energy tables from a file.
    ia_params and force/energy tables are then communicated to each node
    \warning No checking is performed for the file read!!
    @param part_type_a particle type for which the interaction is defined
    @param part_type_b particle type for which the interaction is defined
    @param filename from which file to fetch the data
    @return <ul>
    <li> 0 on success
    <li> 1 on particle type mismatches
    <li> 2 
    </ul>
*/
MDINLINE int tabulated_set_params(int part_type_a, int part_type_b, char* filename)
{
  IA_parameters *data, *data_sym;
  FILE* fp;
  int npoints;
  double minval,minval2, maxval, maxval2;
  int i, newsize;
  int token;
  double dummr;
  token = 0;

  make_particle_type_exist(part_type_a);
  make_particle_type_exist(part_type_b);
    
  data     = get_ia_param(part_type_a, part_type_b);
  data_sym = get_ia_param(part_type_b, part_type_a);

  if (!data || !data_sym)
    return 1;

  if (strlen(filename) > MAXLENGTH_TABFILE_NAME-1 )
    return 2;

  /*Open the file containing force and energy tables */
  fp = fopen( filename , "r");
  if ( !fp )
    return 3;

  /*Look for a line starting with # */
  while ( token != EOF) {
    token = fgetc(fp);
    if ( token == 35 ) { break; } // magic number for # symbol
  }
  if ( token == EOF ) {
    fclose(fp);
    return 4;
  }

  /* First read two important parameters we read in the data later*/
  fscanf( fp , "%d ", &npoints);
  fscanf( fp, "%lf ", &minval);
  fscanf( fp, "%lf ", &maxval);

  // Set the newsize to the same as old size : only changed if a new force table is being added.
  newsize = tabulated_forces.max;

  if ( data->TAB_npoints == 0){
    // A new potential will be added so set the number of points, the startindex and newsize
    data->TAB_npoints    = data_sym->TAB_npoints    = npoints;
    data->TAB_startindex = data_sym->TAB_startindex = tabulated_forces.max;
    newsize += npoints;
  } else {
    // We have existing data for this pair of monomer types check array sizing
    if ( data->TAB_npoints != npoints ) {
      fclose(fp);
      return 5;
    }
  }

  /* Update parameters symmetrically */
  data->TAB_maxval    = data_sym->TAB_maxval    = maxval;
  data->TAB_minval    = data_sym->TAB_minval    = minval;
  strcpy(data->TAB_filename,filename);
  strcpy(data_sym->TAB_filename,filename);

  /* Calculate dependent parameters */
  maxval2 = maxval*maxval;
  minval2 = minval*minval;
  data->TAB_maxval2 = data_sym->TAB_maxval2 = maxval2;
  data->TAB_minval2 = data_sym->TAB_minval2 = minval2;
  data->TAB_stepsize = data_sym->TAB_stepsize = (maxval-minval)/(double)(data->TAB_npoints - 1);


  /* Allocate space for new data */
  realloc_doublelist(&tabulated_forces,newsize);
  realloc_doublelist(&tabulated_energies,newsize);

  /* Read in the new force and energy table data */
  for (i =0 ; i < npoints ; i++)
    {
      fscanf(fp,"%lf",&dummr);
      fscanf(fp,"%lf", &(tabulated_forces.e[i+data->TAB_startindex]));
      fscanf(fp,"%lf", &(tabulated_energies.e[i+data->TAB_startindex]));
    }

  fclose(fp);

  /* broadcast interaction parameters including force and energy tables*/
  mpi_bcast_ia_params(part_type_a, part_type_b);
  mpi_bcast_ia_params(part_type_b, part_type_a);

  if (tab_force_cap != -1.0) {
    mpi_tab_cap_forces(tab_force_cap);}
  return 0;
}

/** Bonded tabulated potentials: Reads tabulated parameters and force
    and energy tables from a file.  ia_params and force/energy tables
    are then communicated to each node \warning No checking is
    performed for the file read!! */
MDINLINE int bonded_tabulated_set_params(int bond_type, int tab_type, char * filename) 
{
  int i, token = 0, size;
  double dummr;
  FILE* fp;

  if(bond_type < 0)
    return 1;
  
  make_bond_type_exist(bond_type);
  
  fp = fopen( filename , "r");
  if ( !fp )
    return 2;
  
  /*Look for a line starting with # */
  while ( token != EOF) {
    token = fgetc(fp);
    if ( token == 35 ) { break; } // magic number for # symbol
  }
  if ( token == EOF ) { 
    fclose(fp);
    return 4;
  }

  /* set types */
  bonded_ia_params[bond_type].type       = BONDED_IA_TABULATED;
  bonded_ia_params[bond_type].p.tab.type = tab_type;

  /* set number of interaction partners */
  if(tab_type == TAB_BOND_LENGTH)   bonded_ia_params[bond_type].num = 1;
  if(tab_type == TAB_BOND_ANGLE)    bonded_ia_params[bond_type].num = 2;
  if(tab_type == TAB_BOND_DIHEDRAL) bonded_ia_params[bond_type].num = 3;

  /* copy filename */
  size = strlen(filename);
  bonded_ia_params[bond_type].p.tab.filename = (char*)malloc((size+1)*sizeof(char));
  strcpy(bonded_ia_params[bond_type].p.tab.filename,filename);

  /* read basic parameters from file */
  fscanf( fp , "%d ", &size);
  bonded_ia_params[bond_type].p.tab.npoints = size;
  fscanf( fp, "%lf ", &bonded_ia_params[bond_type].p.tab.minval);
  fscanf( fp, "%lf ", &bonded_ia_params[bond_type].p.tab.maxval);

  /* Check interval for angle and dihedral potentials.  With adding
     ROUND_ERROR_PREC to the upper boundary we make sure, that during
     the calculation we do not leave the defined table!
  */
  if(tab_type == TAB_BOND_ANGLE ) {
    if( bonded_ia_params[bond_type].p.tab.minval != 0.0 || 
	abs(bonded_ia_params[bond_type].p.tab.maxval-PI) > 1e-5 ) {
      fclose(fp);
      return 4;
    }
    bonded_ia_params[bond_type].p.tab.maxval = PI+ROUND_ERROR_PREC;
  }
  /* check interval for angle and dihedral potentials */
  if(tab_type == TAB_BOND_DIHEDRAL ) {
    if( bonded_ia_params[bond_type].p.tab.minval != 0.0 || 
	abs(bonded_ia_params[bond_type].p.tab.maxval-(2*PI)) > 1e-5 ) {
      fclose(fp);
      return 5;
    }
    bonded_ia_params[bond_type].p.tab.maxval = (2*PI)+ROUND_ERROR_PREC;
  }


				    

  /* calculate dependent parameters */
  bonded_ia_params[bond_type].p.tab.invstepsize = (double)(size-1)/(bonded_ia_params[bond_type].p.tab.maxval-bonded_ia_params[bond_type].p.tab.minval);

  /* allocate force and energy tables */
  bonded_ia_params[bond_type].p.tab.f = (double*)malloc(size*sizeof(double));
  bonded_ia_params[bond_type].p.tab.e = (double*)malloc(size*sizeof(double));

  /* Read in the new force and energy table data */
  for (i =0 ; i < size ; i++) {
      fscanf(fp,"%lf", &dummr);
      fscanf(fp,"%lf", &bonded_ia_params[bond_type].p.tab.f[i]);
      fscanf(fp,"%lf", &bonded_ia_params[bond_type].p.tab.e[i]);
  }
  fclose(fp);

  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the tabulated bonded potential
MDINLINE int inter_parse_bonded_tabulated(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  int tab_type = TAB_UNKNOWN;

  if (argc < 3 ) {
    Tcl_AppendResult(interp, "tabulated needs two string parameter: "
		     "<type> <filename>", (char *) NULL);
    return (TCL_ERROR);
  }  

  if (ARG_IS_S(1,"bond"))     tab_type = TAB_BOND_LENGTH;
  if (ARG_IS_S(1,"angle"))    tab_type = TAB_BOND_ANGLE;
  if (ARG_IS_S(1,"dihedral")) tab_type = TAB_BOND_DIHEDRAL;
  if (tab_type == TAB_UNKNOWN) {
    Tcl_AppendResult(interp, "Unknown type of bonded tabulated interaction. Should be: "
		     "\"bond\" or \"angle\" or \"dihedral\"", (char *) NULL);
    return (TCL_ERROR);
  }

  switch (bonded_tabulated_set_params(bond_type, tab_type, argv[2])) {
  case 1:
    Tcl_AppendResult(interp, "illegal bond type", (char *)NULL);
    return TCL_ERROR;
  case 2:
    Tcl_AppendResult(interp, "cannot open \"", argv[2], "\"", (char *)NULL);
    return TCL_ERROR;
  case 3:
    Tcl_AppendResult(interp, "attempt to read file \"", argv[2], "\" failed."
		     "Could not find start the start token <#>", (char *)NULL);
    return TCL_ERROR;
  case 4:
    Tcl_AppendResult(interp, "bond angle potential has to be defined in the interval 0 to pi", (char *)NULL);
    return TCL_ERROR;
  case 5:
    Tcl_AppendResult(interp, "bond angle potential has to be defined in the interval 0 to 2pi", (char *)NULL);
    return TCL_ERROR;
  default:
    return TCL_OK;
  }
}

/// parser for the force cap
MDINLINE int inter_parse_tabforcecap(Tcl_Interp * interp, int argc, char ** argv)
{
  char buffer[TCL_DOUBLE_SPACE];


  if (argc == 0) {
    if (tab_force_cap == -1.0)
      Tcl_AppendResult(interp, "tabforcecap individual", (char *) NULL);
    else {
      Tcl_PrintDouble(interp, tab_force_cap, buffer);
      Tcl_AppendResult(interp, "tabforcecap ", buffer, (char *) NULL);
    }
    return TCL_OK;
  }

  if (argc > 1) {
    Tcl_AppendResult(interp, "inter tabforcecap takes at most 1 parameter",
		     (char *) NULL);      
    return TCL_ERROR;
  }
  
  if (ARG0_IS_S("individual"))
      tab_force_cap = -1.0;
  else if (! ARG0_IS_D(tab_force_cap) || tab_force_cap < 0) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "force cap must be a nonnegative double value or \"individual\"",
		     (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(tabforcecap_set_params(tab_force_cap),
	      "If you can read this, you should change it. (Use the source Luke!)");
  return TCL_ERROR;
}

MDINLINE int tab_parser(Tcl_Interp * interp,
			int part_type_a, int part_type_b,
			int argc, char ** argv)
{
  char *filename = NULL;

  /* tabulated interactions should supply a file name for a file containing
     both force and energy profiles as well as number of points, max
     values etc.
  */
  if (argc < 2) {
    Tcl_AppendResult(interp, "tabulated potentials require a filename: "
		     "<filename>",
		     (char *) NULL);
    return 0;
  }

  /* copy tabulated parameters */
  filename = argv[1];

  switch (tabulated_set_params(part_type_a, part_type_b, filename)) {
  case 1:
    Tcl_AppendResult(interp, "particle types must be non-negative", (char *) NULL);
    return 0;
  case 2:
    Tcl_AppendResult(interp, "the length of the filename must be less than 256 characters,"
		     "but is \"", filename, "\"", (char *)NULL);
    return 0;
  case 3:
    Tcl_AppendResult(interp, "cannot open \"", filename, "\"", (char *)NULL);
    return 0;
  case 4:
    Tcl_AppendResult(interp, "attempt to read file \"", filename,
		     "\" failed, could not find start the start token <#>", (char *)NULL);
    return 0;
  case 5:
    Tcl_AppendResult(interp, "number of data points does not match the existing table", (char *)NULL);
    return 0;
    
  }
  return 2;
}

/** Add a non-bonded pair force by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void add_tabulated_pair_force(Particle *p1, Particle *p2, IA_parameters *ia_params,
				       double d[3], double dist, double force[3])
{
  double phi, dindex, fac;
  int tablepos, table_start,j;
  double rescaled_force_cap = tab_force_cap/dist;
  double maxval = ia_params->TAB_maxval;
  double minval = ia_params->TAB_minval;

  fac = 0.0;

  if ( maxval > 0 ) {
    if ( dist < maxval){ 
      table_start = ia_params->TAB_startindex;
      dindex = (dist-minval)/ia_params->TAB_stepsize;
      tablepos = (int)(floor(dindex));  

      if ( dist > minval ) {
       phi = dindex - tablepos;	  
       fac = tabulated_forces.e[table_start + tablepos]*(1-phi) + tabulated_forces.e[table_start + tablepos+1]*phi;
      }
      else {
	/* Use an extrapolation beyond the table */
	if ( dist > 0 ) {
	  tablepos = 0;
	  phi = dindex - tablepos;	  
	  fac = (tabulated_forces.e[table_start]*minval)*(1-phi) + 
	    (tabulated_forces.e[table_start+1]*(minval+ia_params->TAB_stepsize))*phi;
	  fac = fac/dist;
	}
	else { /* Particles on top of each other .. leave fac as 0.0 */
	}
      }
      
      if ( rescaled_force_cap < fac && tab_force_cap > 0.0) {
	fac = rescaled_force_cap;
      }
    }
    for(j=0;j<3;j++)
      force[j] += fac * d[j];
  }
}

/** Add a non-bonded pair energy by linear interpolation from a table.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double tabulated_pair_energy(Particle *p1, Particle *p2, IA_parameters *ia_params,
				      double d[3], double dist) {
  double phi, dindex;
  int tablepos, table_start;
  double x0, b;
  

  if ( dist < ia_params->TAB_maxval){ 
    dindex = (dist-ia_params->TAB_minval)/ia_params->TAB_stepsize;
    tablepos = (int)(floor(dindex)); 
    table_start = ia_params->TAB_startindex;

    if ( tablepos < 0 ){
      /* For distances smaller than the tabulated minimim take quadratic
	 extrapolation from first two values of the force table! 
	 This corresponds to the linear extrapolation of the force at that point.
	 This sould not occur too often, since it is quite expensive!
      */
      tablepos = 0;   
      b = (tabulated_forces.e[table_start + tablepos + 1]-tabulated_forces.e[table_start + tablepos])/ia_params->TAB_stepsize;
      x0 = ia_params->TAB_minval-tabulated_forces.e[table_start + tablepos]/b;
      return ( (tabulated_energies.e[table_start + tablepos]
		+ 0.5*b*SQR(ia_params->TAB_minval-x0))
	       - 0.5*b*SQR(dist-x0) );
    }

    phi = (dindex - tablepos);
 
    return  tabulated_energies.e[table_start + tablepos]*(1-phi) 
      + tabulated_energies.e[table_start + tablepos+1]*phi;
  }
  return 0.0;
}


/** check the tabulated forcecap to see that it is sensible \warning
    This routine will probably give strange results if forcecap is
    applied before the table is loaded.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE void check_tab_forcecap(double force_cap)
{
  int i,j,startindex;
  IA_parameters *params;

  for(i=0; i<n_particle_types; i++) {
    for(j=0; j<n_particle_types; j++) {
      params = get_ia_param(i,j);
      startindex = params->TAB_startindex;
      if ( tabulated_forces.max < (params->TAB_npoints + startindex )) { /* Make sure forces are initialized */
	if(force_cap > 0.0 && params->TAB_maxval > 0.0 && tabulated_forces.e[startindex] > force_cap) {
	  for ( i = 0 ; i < params->TAB_npoints ; i++) {
	    if ( tabulated_forces.e[startindex + i] < force_cap ) {
	      return; /* Everything is OK nothing to say :) */
	    }	  
	  }
	  if ( i == params->TAB_npoints - 1) {
	    tab_force_cap = -1.0;
	    /* Force cap is below all known forces .. turn force capping off */
	  }
	}    
	if ( force_cap > tabulated_forces.e[startindex] ) {
	  fprintf(stderr,"calc_tab_cap_radii: Capped region is outside the force table");
	  
	  if ( tabulated_forces.e[startindex] < tabulated_forces.e[startindex+1] ) {
	    fprintf(stderr,"Extrapolation does not make sense outside this table \n");
	    errexit();
	  }
	  
	  fprintf(stderr,", will extrapolate the force outside table \n");
	  fprintf(stderr,"fc: %f \n",tab_force_cap);	
	}
      }
    }
  }
}

/* BONDED INTERACTIONS */


/** Force factor lookup in a force table for bonded interactions (see
    \ref Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.h).*/
MDINLINE double bonded_tab_force_lookup(double val, Bonded_ia_parameters *iaparams)
{
#ifdef TABULATED
  int    ind;
  double dind;
  
  dind = (val - iaparams->p.tab.minval)*iaparams->p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return  iaparams->p.tab.f[ind]*(1.0-dind) + iaparams->p.tab.f[ind+1]*dind;
#else
  return 0.0;
#endif
}

/** Energy lookup in a energy table for bonded interactions (see \ref
    Bonded_ia_parameters). The force is calculated by linear
    interpolation between the closest tabulated values. There is no
    check for the upper bound! 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE double bonded_tab_energy_lookup(double val, Bonded_ia_parameters *iaparams)
{
  int ind;
  double dind;
  
  dind = (val - iaparams->p.tab.minval)*iaparams->p.tab.invstepsize;

  if( dind < 0.0 ) ind = 0;
  else ind = (int)dind;

  dind = dind - ind;
  /* linear interpolation between data points */
  return iaparams->p.tab.e[ind]*(1.0-dind) + iaparams->p.tab.e[ind+1]*dind;
}

/** Calculate a tabulated bond length force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1 and p2 and add it
    to the particle forces. The force acts in the direction of the
    connecting vector between the particles. For distances smaller
    than the tabulated range it uses a linear extrapolation based on
    the first two tabulated force values.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int calc_tab_bond_force(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double force[3]) 
{
  int i;
  double fac, dist = sqrt(sqrlen(dx));

  if(dist > iaparams->p.tab.maxval)
    return 1;

  fac = bonded_tab_force_lookup(dist, iaparams)/dist;
  
  for(i=0;i<3;i++)
    force[i] = -fac*dx[i];

  ONEPART_TRACE(if(p1->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],p2->p.identity,dist,fac));
  ONEPART_TRACE(if(p2->p.identity==check_id) fprintf(stderr,"%d: OPT: TAB BOND f = (%.3e,%.3e,%.3e) with part id=%d at dist %f fac %.3e\n",this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],p1->p.identity,dist,fac));

  return 0;
}

/** Calculate and return a tabulated bond length energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1 and
    p2. For distances smaller than the tabulated range it uses a
    quadratic extrapolation based on the first two tabulated force
    values and the first tabulated energy value. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int tab_bond_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy) 
{
  double dist = sqrt(sqrlen(dx));

  if(dist > iaparams->p.tab.maxval)
    return 1;

  /* For distances smaller than the tabulated minimim take quadratic
     extrapolation from first two values of the force table! 
     This corresponds to the linear extrapolation of the force at that point.
     This sould not occur too often, since it is quite expensive!
  */
  if( dist <  iaparams->p.tab.minval) {
    double x0, b;
    b = (iaparams->p.tab.f[1]-iaparams->p.tab.f[0])*iaparams->p.tab.invstepsize;
    x0 = iaparams->p.tab.minval - iaparams->p.tab.f[0]/b;
    *_energy = ( (iaparams->p.tab.e[0] + 0.5*b*SQR(iaparams->p.tab.minval-x0)) -
		 0.5*b*SQR(dist-x0) );
  }
  else
    *_energy = bonded_tab_energy_lookup(dist, iaparams);

  return 0;
}


/** Calculate a tabulated bond angle force with number type_num (see
    \ref Bonded_ia_parameters) between particles p_left, p_mid and
    p_right and add it to the particle forces. The force on p_left and
    p_right acts perpendicular to the connecting vector between the
    particle and p_mid and in the plane defined by the three
    particles. The force on the middle particle balances the other two
    forces. The forces are scaled with the invers length of the
    connecting vectors. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int calc_tab_angle_force(Particle *p_mid, Particle *p_left, 
				  Particle *p_right, Bonded_ia_parameters *iaparams,
				  double force1[3], double force2[3])
{
    double cosine, phi, invsinphi, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  phi = acos(cosine);
  invsinphi = sin(phi);
  if(invsinphi < TINY_SIN_VALUE) invsinphi = TINY_SIN_VALUE;
  invsinphi = 1.0/invsinphi;
  /* look up force factor */
  fac = bonded_tab_force_lookup(phi, iaparams);
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j])*invsinphi * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j])*invsinphi * d2i;
    force1[j] = (f1-f2);
    force2[j] = -f1;
  }

  return 0;
}

/** Calculate and return tabulated bond angle energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p_left,
    p_mid and p_right. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int tab_angle_energy(Particle *p_mid, Particle *p_left, 
			      Particle *p_right, Bonded_ia_parameters *iaparams,
			      double *_energy)
{
  double phi, vec1[3], vec2[3], vl1, vl2; 

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  vl1 = sqrt(sqrlen(vec1));
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  vl2 = sqrt(sqrlen(vec2));
  /* calculate phi */
  phi = acos( scalar(vec1, vec2) / (vl1*vl2) );

  *_energy = bonded_tab_energy_lookup(phi, iaparams);

  return 0;
}

/** Calculate a tabulated dihedral force with number type_num (see
    \ref Bonded_ia_parameters) between particles p1. p2, p3 and p4 and
    add it to the particle forces. This function is not tested yet.
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int calc_tab_dihedral_force(Particle *p2, Particle *p1,
				     Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
				     double force2[3], double force1[3], double force3[3])
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
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if ( phi == -1.0 ) { 
    for(i=0;i<3;i++) { force1[i] = 0.0; force2[i] = 0.0; force3[i] = 0.0; }
    return 0;
  }

  /* calculate force components (directions) */
  for(i=0;i<3;i++)  {
    f1[i] = (v23Xv34[i] - cosphi*v12Xv23[i])/l_v12Xv23;;
    f4[i] = (v12Xv23[i] - cosphi*v23Xv34[i])/l_v23Xv34;
  }
  vector_product(v23, f1, v23Xf1);
  vector_product(v23, f4, v23Xf4);
  vector_product(v34, f4, v34Xf4);
  vector_product(v12, f1, v12Xf1);

  /* table lookup */
  fac = bonded_tab_force_lookup(phi, iaparams);

  /* store dihedral forces */
  for(i=0;i<3;i++) {
      force1[i] = fac*v23Xf1[i];
      force2[i] = fac*(v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
      force3[i] = fac*(v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }

  return 0;
}

/** Calculate and return a tabulated dihedral energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p1. p2,
    p3 and p4. This function is not tested yet. 
    Needs feature TABULATED compiled in (see \ref config.h). */
MDINLINE int tab_dihedral_energy(Particle *p2, Particle *p1, 
				 Particle *p3, Particle *p4, Bonded_ia_parameters *iaparams,
				 double *_energy)
{
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23, v23Xv34, &l_v23Xv34, &cosphi, &phi);

  *_energy = bonded_tab_energy_lookup(phi, iaparams);

  return 0;
}
#endif

#endif
