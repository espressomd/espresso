// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
/** \file statistics.c
    This is the place for analysation (so far...).
    Implementation of \ref statistics.h "statistics.h"
*/
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include "statistics.h"
#include "statistics_chain.h"
#include "modes.h"
#include "energy.h"
#include "pressure.h"
#include "communication.h"
#include "debug.h"
#include "grid.h"
#include "parser.h"
#include "particle_data.h"
#include "interaction_data.h"

/** Previous particle configurations (needed for offline analysis and correlation analysis in \ref #analyze) */
double **configs = NULL; int n_configs = 0; int n_part_conf = 0;

/****************************************************************************************
 *                                 helper functions
 ****************************************************************************************/

double min_distance2(double pos1[3], double pos2[3])
{
  double diff[3];
  get_mi_vector(diff, pos1, pos2);
  return sqrlen(diff);
}

static int get_reference_point(Tcl_Interp *interp, int *argc, char ***argv,
			       double pos[3], int *pid)
{
  *pid = -1;
  
  if (*argc < 3) {
    Particle ref;
    if (Tcl_GetInt(interp, (*argv)[0], pid) == TCL_ERROR)
      return TCL_ERROR;
    
    if (get_particle_data(*pid, &ref) != TCL_OK) {
      Tcl_AppendResult(interp, "reference particle does not exist", (char *)NULL);
      return TCL_ERROR;
    }
    pos[0] = ref.r.p[0];
    pos[1] = ref.r.p[1];
    pos[2] = ref.r.p[2];
    
    (*argc)--;
    (*argv)++;

    free_particle(&ref);
    return TCL_OK;
  }
  else {
    if (Tcl_GetDouble(interp, (*argv)[0], &pos[0]) == TCL_ERROR ||
	Tcl_GetDouble(interp, (*argv)[1], &pos[1]) == TCL_ERROR ||
	Tcl_GetDouble(interp, (*argv)[2], &pos[2]) == TCL_ERROR)
      return TCL_ERROR;

    (*argc) -= 3;
    (*argv) += 3;

    return TCL_OK;
  }
  return TCL_ERROR;
}

/****************************************************************************************
 *                                 basic observables calculation
 ****************************************************************************************/

double mindist(IntList *set1, IntList *set2)
{
  double mindist, pt[3];
  int i, j;

  mindist = SQR(box_l[0] + box_l[1] + box_l[2]);

  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles-1; j++) {
    pt[0] = partCfg[j].r.p[0];
    pt[1] = partCfg[j].r.p[1];
    pt[2] = partCfg[j].r.p[2];
    if (!set1 || intlist_contains(set1, partCfg[j].p.type)) {
      for (i=j+1; i<n_total_particles; i++)
	if (!set2 || intlist_contains(set2, partCfg[i].p.type)) {
	  mindist = dmin(mindist, min_distance2(pt, partCfg[i].r.p));
	}
    }
  }
  mindist = sqrt(mindist);
  return mindist;
}

void centermass(int type, double *com)
{
  int i, j,count;

  com[0]=com[1]=com[2]=0.;
  count=0;
   	
  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles; j++) {
    if (type == partCfg[j].p.type) {
      count ++;
      for (i=0; i<3; i++) {
      	com[i] += partCfg[j].r.p[i];
      }
    }
  }
  
  for (i=0; i<3; i++) {
    com[i] /= count;
  }
  return;
}

void gyrationtensor(int type, double *gyrtensor)
{
  int i,j,count;
  double p1[3],com[3];

  count=0;
  updatePartCfg(WITHOUT_BONDS);
  for(i=0;i<9;i++) gyrtensor[i]=0.;
	centermass(type, com);
  for (j=0; j<n_total_particles; j++) {
    if (type == partCfg[j].p.type) {
      count ++;
      for (i=0; i<3; i++) {
      	p1[i] = partCfg[j].r.p[i] - com[i];
      }
    	gyrtensor[0] += p1[0] * p1[0]; 
    	gyrtensor[4] += p1[1] * p1[1];
    	gyrtensor[8] += p1[2] * p1[2];
    	gyrtensor[1] += p1[0] * p1[1];
    	gyrtensor[2] += p1[0] * p1[2]; 
    	gyrtensor[5] += p1[1] * p1[2];
		}
	}
  /* use symmetry */
  gyrtensor[3] = gyrtensor[1]; 
  gyrtensor[6] = gyrtensor[2]; 
  gyrtensor[7] = gyrtensor[5];

  return;
}

void nbhood(double pt[3], double r, IntList *il)
{
  double d[3];
  int i;

  init_intlist(il);
 
  updatePartCfg(WITHOUT_BONDS);

  for (i = 0; i<n_total_particles; i++) {
    get_mi_vector(d, pt, partCfg[i].r.p);
    if (sqrt(sqrlen(d)) < r) {
      realloc_intlist(il, il->n + 1);
      il->e[il->n] = partCfg[i].p.identity;
      il->n++;
    }
  }
}

double distto(double p[3], int pid)
{
  int i;
  double d[3];
  double mindist;

  /* larger than possible */
  mindist=SQR(box_l[0] + box_l[1] + box_l[2]);
  for (i=0; i<n_total_particles; i++) {
    if (pid != partCfg[i].p.identity) {
      get_mi_vector(d, p, partCfg[i].r.p);
      mindist = dmin(mindist, sqrlen(d));
    }
  }
  return sqrt(mindist);
}

void calc_part_distribution(int *p1_types, int n_p1, int *p2_types, int n_p2, 
			    double r_min, double r_max, int r_bins, int log_flag, 
			    double *low, double *dist)
{
  int i,j,t1,t2,ind,cnt=0;
  double inv_bin_width=0.0;
  double min_dist,min_dist2=0.0,start_dist2,act_dist2;

  start_dist2 = SQR(box_l[0] + box_l[1] + box_l[2]);
  /* bin preparation */
  *low = 0.0;
  for(i=0;i<r_bins;i++) dist[i] = 0.0;
  if(log_flag == 1) inv_bin_width = (double)r_bins/(log(r_max)-log(r_min));
  else              inv_bin_width = (double)r_bins / (r_max-r_min);

  /* particle loop: p1_types*/
  for(i=0; i<n_total_particles; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	min_dist2 = start_dist2;
	/* particle loop: p2_types*/
	for(j=0; j<n_total_particles; j++) {
	  if(j != i) {
	    for(t2=0; t2<n_p2; t2++) {
	      if(partCfg[j].p.type == p2_types[t2]) {
		act_dist2 =  min_distance2(partCfg[i].r.p, partCfg[j].r.p);
		if(act_dist2 < min_dist2) { min_dist2 = act_dist2; }
	      }
	    }
	  }
	}
	min_dist = sqrt(min_dist2);
	if(min_dist <= r_max) {
	  if(min_dist >= r_min) {
	    /* calculate bin index */
	    if(log_flag == 1) ind = (int) ((log(min_dist) - log(r_min))*inv_bin_width);
	    else              ind = (int) ((min_dist - r_min)*inv_bin_width);
	    if(ind >= 0 && ind < r_bins) {
	      dist[ind] += 1.0;
	    }
	  }
	  else {
	    *low += 1.0;
	  }
	}
	cnt++;    
      }
    }
  }
  
  /* normalization */
  *low /= (double)cnt;
  for(i=0;i<r_bins;i++) dist[i] /= (double)cnt;
}

void calc_rdf(int *p1_types, int n_p1, int *p2_types, int n_p2, 
	      double r_min, double r_max, int r_bins, double *rdf)
{
  int i,j,t1,t2,ind,cnt=0;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;

  if(n_p1 == n_p2) {
    for(i=0;i<n_p1;i++) 
      if( p1_types[i] != p2_types[i] ) mixed_flag=1;
  }
  else mixed_flag=1;

  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  for(i=0;i<r_bins;i++) rdf[i] = 0.0;
  /* particle loop: p1_types*/
  for(i=0; i<n_total_particles; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	/* distinguish mixed and identical rdf's */
	if(mixed_flag == 1) start = 0;
	else                start = (i+1);
	/* particle loop: p2_types*/
	for(j=start; j<n_total_particles; j++) {
	  for(t2=0; t2<n_p2; t2++) {
	    if(partCfg[j].p.type == p2_types[t2]) {
	      dist = min_distance(partCfg[i].r.p, partCfg[j].r.p);
	      if(dist > r_min && dist < r_max) {
		ind = (int) ( (dist - r_min)*inv_bin_width );
		rdf[ind]++;
	      }
	      cnt++;
	    }
	  }
	}
      }
    }
  }

  /* normalization */
  volume = box_l[0]*box_l[1]*box_l[2];
  for(i=0; i<r_bins; i++) {
    r_in       = i*bin_width + r_min; 
    r_out      = r_in + bin_width;
    bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
    rdf[i] *= volume / (bin_volume * cnt);
  }
}

/****************************************************************************************
 *                                 config storage functions
 ****************************************************************************************/

void analyze_append() {
  int i;
  n_part_conf = n_total_particles;
  configs = realloc(configs,(n_configs+1)*sizeof(double *));
  configs[n_configs] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs][3*i]   = partCfg[i].r.p[0];
    configs[n_configs][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs][3*i+2] = partCfg[i].r.p[2];
  }
  n_configs++;
}

void analyze_push() {
  int i;
  n_part_conf = n_total_particles;
  free(configs[0]);
  for(i=0; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  configs[n_configs-1] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs-1][3*i]   = partCfg[i].r.p[0];
    configs[n_configs-1][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs-1][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_replace(int ind) {
  int i;
  n_part_conf = n_total_particles;
  for(i=0; i<n_part_conf; i++) {
    configs[ind][3*i]   = partCfg[i].r.p[0];
    configs[ind][3*i+1] = partCfg[i].r.p[1];
    configs[ind][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_remove(int ind) {
  int i;
  free(configs[ind]);
  for(i=ind; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  n_configs--;
  configs = realloc(configs,n_configs*sizeof(double *));
  if (n_configs == 0) n_part_conf = 0;
}

void analyze_configs(double *tmp_config, int count) {
  int i;
  n_part_conf = count;
  configs = realloc(configs,(n_configs+1)*sizeof(double *));
  configs[n_configs] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs][3*i]   = tmp_config[3*i];
    configs[n_configs][3*i+1] = tmp_config[3*i+1];
    configs[n_configs][3*i+2] = tmp_config[3*i+2];
  }
  n_configs++;
}

/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded, int n_non_bonded,
			       int n_coulomb, int c_size)
{
  int i, total = c_size*(n_pre + n_bonded_ia + n_non_bonded + n_coulomb);

  realloc_doublelist(&(stat->data), stat->data.n = total);
  stat->chunk_size = c_size;
  stat->n_coulomb    = n_coulomb;
  stat->n_non_bonded = n_non_bonded;
  stat->bonded     = stat->data.e + c_size*n_pre;
  stat->non_bonded = stat->bonded + c_size*n_bonded_ia;
  stat->coulomb    = stat->non_bonded + c_size*n_non_bonded;

  for(i = 0; i < total; i++)
    stat->data.e[i] = 0.0;
}

void invalidate_obs()
{ 
  total_energy.init_status = 0;
  total_pressure.init_status = 0;
}

/****************************************************************************************
 *                                 basic observables parsing
 ****************************************************************************************/

static int parse_modes2d(Tcl_Interp *interp, int argc, char **argv)
{
  STAT_TRACE(fprintf(stderr,"%d,parsing modes2d \n",this_node);)
  /* 'analyze modes2d [setgrid <xdim> <ydim> <zdim>] [setstray <stray_cut_off>]]' */
  char buffer[TCL_DOUBLE_SPACE];
  int i,j,change ;
  fftw_complex* result;

  change = 0;


  if (n_total_particles <= 2) {
    Tcl_AppendResult(interp, "(not enough particles for mode analysis)",
		     (char *)NULL);
    return (TCL_OK);
  }

  while (argc > 0)
    {
      if ( ARG0_IS_S("setgrid") ) { 
	if ( !ARG_IS_I(1,mode_grid_3d[0]) || !ARG_IS_I(2,mode_grid_3d[1]) || !ARG_IS_I(3,mode_grid_3d[2]) ) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp,"usage: analyze modes2d [setgrid <xdim> <ydim> <zdim>] [setstray <stray_cut_off>]", (char *)NULL);
	  return (TCL_ERROR);
	}
	STAT_TRACE(fprintf(stderr,"%d,setgrid has args %d,%d,%d \n",this_node,mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]));
	change = 4;
	/* Update global parameters */
	map_to_2dgrid();
	mode_grid_changed = 1;
	
      }
      if ( ARG0_IS_S("setstray") ) { 
	if ( !ARG_IS_D(1,stray_cut_off) ) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp,"usage: analyze modes2d [setgrid <xdim> <ydim> <zdim>] [setstray <stray_cut_off>]", (char *)NULL);
	  return (TCL_ERROR);
	}
	STAT_TRACE(fprintf(stderr,"%d,setgrid has args %d,%d,%d \n",this_node,mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]));
	change = 2;
      }
      argc -= change;
      argv += change;
      STAT_TRACE(fprintf(stderr,"%d,argc = %d \n",this_node, argc);)

	
    }


  result = malloc((mode_grid_3d[ydir]/2+1)*(mode_grid_3d[xdir])*sizeof(fftw_complex));

  if (!modes2d(result)) {
    fprintf(stderr,"%d,mode analysis failed \n",this_node);
  }
  else {    STAT_TRACE(fprintf(stderr,"%d,mode analysis done \n",this_node));}
  

  Tcl_AppendResult(interp, "{ Modes } { ", (char *)NULL);
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    Tcl_AppendResult(interp, " { ", (char *)NULL);
    for ( j = 0 ; j < mode_grid_3d[ydir]/2 + 1 ; j++) {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      Tcl_PrintDouble(interp,result[j+i*(mode_grid_3d[ydir]/2+1)].re,buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
      Tcl_PrintDouble(interp,result[j+i*(mode_grid_3d[ydir]/2+1)].im,buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " } ", (char *)NULL);
      }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
  }


  Tcl_AppendResult(interp, " } ", (char *)NULL);

  free(result);

  return TCL_OK;

}

static int parse_mindist(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze mindist [<type_list_a> <type_list_b>]' */
  double result;
  char buffer[TCL_DOUBLE_SPACE];
  IntList p1, p2;

  init_intlist(&p1); init_intlist(&p2);

  if (n_total_particles <= 1) {
    Tcl_AppendResult(interp, "(not enough particles)",
		     (char *)NULL);
    return (TCL_OK);
  }
  if (argc == 0)
    result = mindist(NULL, NULL);
  else {
    /* parse arguments */
    if (argc < 2) {
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }

    if (!ARG0_IS_INTLIST(p1)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }
    if (!ARG1_IS_INTLIST(p2)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }
    result = mindist(&p1, &p2);
  }

  Tcl_PrintDouble(interp, result, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}
static int parse_centermass(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze centermass [<type>]' */
  double com[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1;
  
  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze centermass [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze centermass [<type>]", (char *)NULL);
      return (TCL_ERROR);
  }
  
  centermass(p1, com);
  
  sprintf(buffer,"%f %f %f",com[0],com[1],com[2]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_gyrationtensor(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze gyrationtensor [<type>]' */
  double gyrtensor[9];
  char buffer[9*TCL_DOUBLE_SPACE+9];
  int p1;

  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze gyrationtensor [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze gyrationtensor [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  gyrationtensor(p1, gyrtensor);
  
  sprintf(buffer,"%f %f %f %f %f %f %f %f %f",
    gyrtensor[0],gyrtensor[1],gyrtensor[2],gyrtensor[3],gyrtensor[4],
    gyrtensor[5],gyrtensor[6],gyrtensor[7],gyrtensor[8]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_find_principal_axis(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze find_principal_axis [<type0>]' */
  double gyrtensor[9],eva[3],eve[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1,i,j;

  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze find_principal_axis [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze find_principal_axis [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  gyrationtensor(p1, gyrtensor);
  i=calc_eigenvalues_3x3(gyrtensor, eva);
  
  sprintf(buffer,"{eigenval eigenvector} ");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  for (j= 0; j < 3; j++) {
  	i=calc_eigenvector_3x3(gyrtensor,eva[j],eve);
    sprintf(buffer," { %f { %f %f %f } }",eva[j],eve[0],eve[1],eve[2]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  return TCL_OK;
}

static int parse_nbhood(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze nbhood { <partid> | <posx> <posy> <posz> } <r_catch>' */
  int p, i;
  double pos[3];
  double r_catch;
  char buffer[TCL_INTEGER_SPACE + 2];  
  IntList il;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  get_reference_point(interp, &argc, &argv, pos, &p);

  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: nbhood { <partid> | <posx> <posy> <posz> } <r_catch>",
		     (char *)NULL);
    return TCL_ERROR;
  }

  if (!ARG0_IS_D(r_catch))
    return (TCL_ERROR);

  updatePartCfg(WITHOUT_BONDS);

  nbhood(pos, r_catch, &il);
  
  for (i = 0; i < il.n; i++) {
    sprintf(buffer, "%d ", il.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  realloc_intlist(&il, 0);
  return (TCL_OK);
}

static int parse_distto(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze distto { <part_id> | <posx> <posy> <posz> }' */
  double result;
  int p;
  double pos[3];
  char buffer[TCL_DOUBLE_SPACE];  

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  get_reference_point(interp, &argc, &argv, pos, &p);
  if (argc != 0) {
    Tcl_AppendResult(interp, "usage: distto { <partid> | <posx> <posy> <posz> }",
		     (char *)NULL);
    return TCL_ERROR;
  }

  updatePartCfg(WITHOUT_BONDS);

  result = distto(pos, p);

  Tcl_PrintDouble(interp, result, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}


static int parse_distribution(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze distribution { <part_type_list_a> } { <part_type_list_b> } [<r_min> [<r_max> [<r_bins> [<log_flag> [<int_flag>]]]]]' */
  /*********************************************************************************************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+3*TCL_INTEGER_SPACE+256];
  IntList p1,p2;
  double r_min=0, r_max=-1.0;
  int r_bins=0, log_flag=0, int_flag=0;
  int i;
  double *distribution, low;

  init_intlist(&p1); init_intlist(&p2);

  if (argc < 2) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_INTLIST(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG1_IS_INTLIST(p2)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }

  argc -= 2; argv += 2;

  if( argc>0 ) { if (!ARG0_IS_D(r_min)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_D(r_max)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(r_bins)) return (TCL_ERROR);   argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(log_flag)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(int_flag)) return (TCL_ERROR); argc--; argv++; }

  /* if not given set defaults */
  if(r_max == -1.) r_max = min_box_l/2.0;
  if(r_bins < 0 )  r_bins = n_total_particles / 20;

  /* give back what you do */
  Tcl_AppendResult(interp, "{ analyze distribution { ", (char *)NULL);
  for(i=0; i<p1.max; i++) { 
    sprintf(buffer,"%d ",p1.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp,"} { ", (char *)NULL);
  for(i=0; i<p2.max; i++) {
    sprintf(buffer,"%d ",p2.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer,"} %f %f %d %d %d",r_min,r_max,r_bins,log_flag,int_flag);
  Tcl_AppendResult(interp, buffer," }", (char *)NULL);
  /* some sanity checks */
  if(r_min < 0.0 || (log_flag==1 && r_min ==0.0 )) return TCL_ERROR;
  if(r_max <= r_min) return TCL_ERROR;
  if(r_bins < 1) return TCL_ERROR;
  /* calculate distribution */
  distribution = malloc(r_bins*sizeof(double));
  updatePartCfg(WITHOUT_BONDS);
  calc_part_distribution(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, log_flag,&low,distribution);
  if(int_flag==1) {
    distribution[0] += low;
    for(i=0; i<r_bins-1; i++) distribution[i+1] += distribution[i]; 
  }
  /* append result */
  {
    double log_fac=0.0, bin_width=0.0, r=0.0;
    if(log_flag == 1) {
      log_fac       = pow((r_max/r_min),(1.0/(double)r_bins));
      r = r_min * sqrt(log_fac);
    } 
    else {
      bin_width     = (r_max-r_min) / (double)r_bins;
      r = r_min + bin_width/2.0;
    }
    Tcl_AppendResult(interp, " {\n", (char *)NULL);
    for(i=0; i<r_bins; i++) {
      sprintf(buffer,"%f %f",r,distribution[i]);
      Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      if(log_flag == 1) r *= log_fac; else r += bin_width;
    }
    Tcl_AppendResult(interp, "}\n", (char *)NULL);
  }
  free(distribution);
  return (TCL_OK);
}

static int parse_rdf(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze rdf' (radial distribution function) */
  /************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE+256];
  IntList p1,p2;
  double r_min=0, r_max=-1.0;
  int r_bins=-1, i;
  double *rdf;

  init_intlist(&p1); init_intlist(&p2);

  if (argc < 2) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze rdf [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_INTLIST(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze rdf [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG1_IS_INTLIST(p2)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze rdf [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }
  argc-=2; argv+=2;

  if( argc>0 ) { if (!ARG0_IS_D(r_min)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_D(r_max)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(r_bins)) return (TCL_ERROR); argc--; argv++; }

  /* if not given use default */
  if(r_max == -1.0) r_max = min_box_l/2.0;
  if(r_bins == -1) r_bins = n_total_particles / 20;

  /* give back what you do */
  Tcl_AppendResult(interp, "{ analyze rdf { ", (char *)NULL);
  for(i=0; i<p1.max; i++) { 
    sprintf(buffer,"%d ",p1.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp,"} { ", (char *)NULL);
  for(i=0; i<p2.max; i++) {
    sprintf(buffer,"%d ",p2.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer,"} %f %f %d",r_min,r_max,r_bins);
  Tcl_AppendResult(interp, buffer," }", (char *)NULL);
  rdf = malloc(r_bins*sizeof(double));
  updatePartCfg(WITHOUT_BONDS);
  calc_rdf(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, rdf);
  /* append result */
  {
    double bin_width=0.0, r=0.0;
    bin_width     = (r_max-r_min) / (double)r_bins;
    r = r_min + bin_width/2.0;
    Tcl_AppendResult(interp, " {\n", (char *)NULL);
    for(i=0; i<r_bins; i++) {
      sprintf(buffer,"%f %f",r,rdf[i]);
      Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      r += bin_width;
    }
    Tcl_AppendResult(interp, "}\n", (char *)NULL);
  }
  free(rdf);
  return (TCL_OK);
}

/****************************************************************************************
 *                                 parser for config storage stuff
 ****************************************************************************************/

static int parse_append(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze append' */
  /********************/
  char buffer[2*TCL_INTEGER_SPACE+256];

  if (argc != 0) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze append", (char *)NULL); return TCL_ERROR; }
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  analyze_append();
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_push(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze push [<size>]' */
  /*****************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i, j;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (argc == 1) { 
    if(!ARG0_IS_I(i)) return (TCL_ERROR);
    if (n_configs < i) analyze_append(); else analyze_push();
    if (n_configs > i) for(j=0; j < n_configs-i; j++) analyze_remove(0);
  }
  else if (argc != 0) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze push [<size>]", (char *)NULL); return TCL_ERROR; }
  else if (n_configs > 0) analyze_push();
  else analyze_append();
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_replace(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze replace <index>' */
  /*****************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i;

  if (argc != 1) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze replace <index>", (char *)NULL); return TCL_ERROR; }
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (!ARG0_IS_I(i)) return (TCL_ERROR);
  if((n_configs == 0) && (i==0)) analyze_append();
  else if ((n_configs == 0) && (i!=0)) {
    Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be replaced!", (char *)NULL); return TCL_ERROR; }
  else if((i < 0) || (i > n_configs-1)) {
    sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
    Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
  else analyze_replace(i);
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_remove(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze remove [<index>]' */
  /******************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i;

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (argc == 0) { for (i = n_configs-1; i >= 0; i--) analyze_remove(i); }
  else if (argc == 1) {
    if (!ARG0_IS_I(i)) return (TCL_ERROR);
    if(n_configs == 0) {
      Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be removed!", (char *)NULL); return TCL_ERROR; }
    else if((i < 0) || (i > n_configs-1)) {
      sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
    analyze_remove(i);
  }
  else {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze remove [<index>]", (char *)NULL); return TCL_ERROR; 
  }
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_configs(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze configs [ { <which> | <configuration> } ]' */
  /*******************************************************/
  char buffer[3*TCL_DOUBLE_SPACE+4*TCL_INTEGER_SPACE+256];
  double *tmp_config;
  int i, j;

  if (argc == 0) {
    for(i=0; i < n_configs; i++) {
      Tcl_AppendResult(interp,"{ ", (char *)NULL);
      for(j=0; j < n_part_conf; j++) {
	sprintf(buffer,"%f %f %f ",configs[i][3*j],configs[i][3*j+1],configs[i][3*j+2]);
	Tcl_AppendResult(interp, buffer,(char *)NULL);
      }
      Tcl_AppendResult(interp,"} ",(char *)NULL);
    }
    return (TCL_OK); }
  else if (argc == 1) {
    if (!ARG0_IS_I(i)) return (TCL_ERROR);
    if ((i<0) || (i>n_configs-1)) {
      sprintf(buffer,"The configs[%d] you requested does not exist, argument must be in [0,%d]!",i,n_configs-1);
      Tcl_AppendResult(interp,buffer,(char *)NULL); return TCL_ERROR; }
    for(j=0; j < n_part_conf; j++) {
      sprintf(buffer,"%f %f %f ",configs[i][3*j],configs[i][3*j+1],configs[i][3*j+2]);
      Tcl_AppendResult(interp, buffer,(char *)NULL);
    }
    return (TCL_OK); }
  else if ((argc == 3*n_part_conf) || (n_part_conf == 0)) {
    if ((n_part_conf == 0) && (argc % 3 == 0)) n_part_conf = argc/3;
    else if (argc != 3*n_part_conf) {
      sprintf(buffer,"Wrong # of args(%d)! Usage: analyze configs [x0 y0 z0 ... x%d y%d z%d]",argc,n_part_conf,n_part_conf,n_part_conf);
      Tcl_AppendResult(interp,buffer,(char *)NULL); return TCL_ERROR; }
    tmp_config = malloc(3*n_part_conf*sizeof(double));
    for(j=0; j < argc; j++)
      if (!ARG_IS_D(j, tmp_config[j])) return (TCL_ERROR);
    analyze_configs(tmp_config, n_part_conf); free(tmp_config);
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK; }
  else { 
    sprintf(buffer,"Wrong # of args(%d)! Usage: analyze configs [x0 y0 z0 ... x%d y%d z%d]",argc,n_part_conf,n_part_conf,n_part_conf);
    Tcl_AppendResult(interp,buffer,(char *)NULL); return TCL_ERROR;
  }
  /* keep the compiler happy */
  return TCL_ERROR;
}

/****************************************************************************************
 *                                 main parser for analyze
 ****************************************************************************************/

int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> ...", (char *)NULL);
    return (TCL_ERROR);
  }

  if (ARG1_IS_S("set"))
    return parse_analyze_set_topology(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("modes2d"))
    return parse_modes2d(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("mindist"))
    return parse_mindist(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("centermass"))
    return parse_centermass(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("gyrationtensor"))
    return parse_gyrationtensor(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("find_principal_axis"))
    return parse_find_principal_axis(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("nbhood"))
    return parse_nbhood(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("distto"))
    return parse_distto(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("energy"))
    return parse_and_print_energy(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("pressure"))
    return parse_and_print_pressure(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("bins"))
    return parse_bins(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("p_IK1"))
    return parse_and_print_p_IK1(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("re"))
    return parse_re(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<re>"))
    return parse_re(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("rg"))
    return parse_rg(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<rg>"))
    return parse_rg(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("rh"))
    return parse_rh(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<rh>"))
    return parse_rh(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("internal_dist"))
    return parse_intdist(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<internal_dist>"))
    return parse_intdist(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("bond_l"))
    return parse_bond_l(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<bond_l>"))
    return parse_bond_l(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("bond_dist"))
    return parse_bond_dist(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<bond_dist>"))
    return parse_bond_dist(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("g123"))
    return parse_g123(interp, 1, argc - 2, argv + 2);    
  else if (ARG1_IS_S("<g1>"))
    return parse_g_av(interp, 1, argc - 2, argv + 2);    
  else if (ARG1_IS_S("<g2>"))
    return parse_g_av(interp, 2, argc - 2, argv + 2);    
  else if (ARG1_IS_S("<g3>"))
    return parse_g_av(interp, 3, argc - 2, argv + 2);
  else if (ARG1_IS_S("formfactor"))
    return parse_formfactor(interp, 0, argc - 2, argv + 2);
  else if (ARG1_IS_S("<formfactor>"))
    return parse_formfactor(interp, 1, argc - 2, argv + 2);
  else if (ARG1_IS_S("distribution"))
    return parse_distribution(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("rdf"))
    return parse_rdf(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("append"))
    return parse_append(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("push"))
    return parse_push(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("replace"))
    return parse_replace(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("remove"))
    return parse_remove(interp, argc - 2, argv + 2);
  else if (ARG1_IS_S("stored")) {
    /* 'analyze stored' */
    /********************/
    char buffer[TCL_INTEGER_SPACE];
    if (argc != 2) {
      Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze stored", (char *)NULL); return TCL_ERROR; 
    }
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  else if (ARG1_IS_S("configs"))
    return parse_configs(interp, argc - 2, argv + 2);
  else {
    /* the default */
    /***************/
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"", argv[1],
		     "\" you requested is not implemented.", (char *)NULL);
    return (TCL_ERROR);
  }
  return (TCL_ERROR);
}
