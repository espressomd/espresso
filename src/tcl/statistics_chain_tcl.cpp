/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file statistics_chain.cpp
    Implementation of \ref statistics_chain.hpp "statistics_chain.h".
*/
#include "statistics.hpp"
#include "parser.hpp"
#include "statistics_chain.hpp"
#include "topology_tcl.hpp"

/****************************************************************************************
 *                                 chain structure commands parsing
 ****************************************************************************************/

int tclcommand_analyze_set_parse_chain_topology(Tcl_Interp *interp, int argc, char **argv)
{
  /* parses a chain topology (e.g. in 'analyze ( rg | <rg> ) [chain start n chains chain length]' , or
     in 'analyze set chains <chain_start> <n_chains> <chain_length>') */
  int m, i, pc;
  
  if (argc < 3) {
    Tcl_AppendResult(interp, "chain structure info consists of <start> <n> <length>", (char *)NULL);    
    return TCL_ERROR;
  }

  if (! (ARG0_IS_I(chain_start) && ARG1_IS_I(chain_n_chains) && ARG_IS_I(2, chain_length)))
    return TCL_ERROR;

  realloc_topology(chain_n_chains);
  pc = 0;
  for (m = 0; m < n_molecules; m++) {
    topology[m].type = 0;
    realloc_intlist(&topology[m].part, topology[m].part.n = chain_length);
    for (i = 0; i < chain_length; i++)
      topology[m].part.e[i] = pc++;
  }
 
  return TCL_OK;
  return tclcommand_analyze_set_parse_topo_part_sync(interp); 
}




/** this function scans the arguments for a description of the chain structure,
    i.e. start of chains, number of chains and chain length. Since this structure
    requires the particles to be sorted, this is performed, too. */

static int tclcommand_analyze_set_parse_chain_topology_check(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc > 0)
    if (tclcommand_analyze_set_parse_chain_topology(interp, argc, argv) != TCL_OK)
      return TCL_ERROR;
  
  if (!sortPartCfg()) {
    Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  return TCL_OK;
}

int tclcommand_analyze_parse_re(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { re | <re> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[4*TCL_DOUBLE_SPACE+4];
  double *re;

  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }

  if (!average)
    calc_re(&re);
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze re' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    calc_re_av(&re);
  }

  sprintf(buffer,"%f %f %f %f",re[0],re[1],re[2],re[3]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(re);
  return (TCL_OK);
}


int tclcommand_analyze_parse_rg(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { rg | <rg> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[4*TCL_DOUBLE_SPACE+4];
  double *rg;
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }
  if (!average)
    calc_rg(&rg);
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rg' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    calc_rg_av(&rg);
  }

  sprintf(buffer,"%f %f %f %f",rg[0],rg[1],rg[2],rg[3]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(rg);
  return (TCL_OK);
}

int tclcommand_analyze_parse_rh(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { rh | <rh> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[2*TCL_DOUBLE_SPACE+2];
  double *rh;
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }
  if (!average)
    calc_rh(&rh);
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rh' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    else
      calc_rh_av(&rh);
  }

  sprintf(buffer,"%f %f",rh[0],rh[1]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(rh);
  return (TCL_OK);
}

int tclcommand_analyze_parse_internal_dist(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { internal_dist | <internal_dist> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[TCL_DOUBLE_SPACE+2];
  int i;
  double *idf;

  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (!average)
    calc_internal_dist(&idf); 
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze internal_dist' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    else
      calc_internal_dist_av(&idf);
  }

  for (i=0; i<chain_length; i++) { 
    sprintf(buffer,"%f ",idf[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }

  free(idf);
  return (TCL_OK);
}

int tclcommand_analyze_parse_bond_l(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { bond_l | <bond_l> } [<chain_start> <n_chains> <chain_length>]' */
  /*****************************************************************************/
  char buffer[4*TCL_DOUBLE_SPACE+4];
  double *bond_l;

  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (!average) calc_bond_l(&bond_l); 
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze bond_l' to only look at current state!", (char *)NULL);
    return TCL_ERROR; }
  else calc_bond_l_av(&bond_l);
  sprintf(buffer,"%f %f %f %f",bond_l[0],bond_l[1],bond_l[2],bond_l[3]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  free(bond_l);
  return (TCL_OK);
}

int tclcommand_analyze_parse_bond_dist(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { bond_dist | <bond_dist> } [index <index>] [<chain_start> <n_chains> <chain_length>]' */
  /***************************************************************************************************/
  char buffer[100 + TCL_DOUBLE_SPACE + 3*TCL_INTEGER_SPACE];
  double *bdf; int ind_n=0, i;

  if (argc >= 2 && ARG_IS_S(0,"index") )
  { 
    if( !ARG_IS_I(1,ind_n) )
    { 
      sprintf(buffer, "\nWrong type for argument");
      Tcl_AppendResult(interp, buffer, (char *)NULL);  
      return (TCL_ERROR); 
    }

    argc-=2; 
    argv+=2; 
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (ind_n < 0 || ind_n > chain_length-1) { 
    sprintf(buffer,"%d!",chain_length-1);
    Tcl_AppendResult(interp, "ERROR: <index> must be between 0 and ", buffer, (char *)NULL);  return TCL_ERROR; }
  if (ind_n >= chain_length/2) ind_n = (chain_length-1) - ind_n;
  if (!average) calc_bond_dist(&bdf,ind_n); 
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze internal_dist' to only look at current state!", (char *)NULL);
    return TCL_ERROR; }
  else calc_bond_dist_av(&bdf,ind_n);
  for (i=0; i<chain_length-ind_n; i++) { 
    sprintf(buffer,"%f ",bdf[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }
  free(bdf);
  return (TCL_OK);
}

	   
int tclcommand_analyze_parse_g123(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze g123 [-init] [<chain_start> <n_chains> <chain_length>]' */
  /********************************************************************/
  char buffer[3*TCL_DOUBLE_SPACE+7];
  int init = 0;
  double g1, g2, g3;

  if (argc > 0 && ARG0_IS_S("-init")) {
    init = 1; argc--; argv++; 
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  
  if (init) { init_g123(); return TCL_OK; }
  if (partCoord_g == NULL || partCM_g == NULL) {
    Tcl_AppendResult(interp, "please call with -init first", (char *)NULL); return TCL_ERROR; }
  if (chain_n_chains != n_chains_g || n_part != n_part_g) {
    fprintf(stderr, "%d %d %d %d\n", chain_n_chains, n_chains_g, n_part, n_part_g);
    Tcl_AppendResult(interp, "initial config has different topology", (char *)NULL);
    return TCL_ERROR;      
  }
  calc_g123(&g1, &g2, &g3);
  sprintf(buffer,"{ %f %f %f }",g1, g2, g3);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}

int tclcommand_analyze_parse_g_av(Tcl_Interp *interp, int what, int argc, char **argv)
{
  /* 'analyze { <g1> | <g2> | <g3> } [<chain_start> <n_chains> <chain_length>]' */
  /******************************************************************************/
  int i, window = 1;
  char buffer[TCL_DOUBLE_SPACE+2];
  double *gx = NULL;
  double weights[3]; weights[0] = weights[1] = weights[2] = 1;

  while (argc >= 1) {
    if (ARG0_IS_S("-sliding")) {
      window = 0;
      argc--; argv++;
    }
    else if (ARG0_IS_S("-weights")) {
      if (argc < 4 || !ARG_IS_D(1,weights[0]) || !ARG_IS_D(2,weights[1]) || !ARG_IS_D(3,weights[2]) )
	return (TCL_ERROR);
      argc-=4; argv+=4;
    }
    else break;
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info or -sliding allowed", (char *)NULL); return TCL_ERROR; }
  if (n_configs == 0) { Tcl_AppendResult(interp, "no configurations found! Use 'analyze append' to save some!", (char *)NULL); return TCL_ERROR; }
  switch (what) {
  case 1:
    calc_g1_av(&gx, window, weights); break;
  case 2:
    calc_g2_av(&gx, window, weights); break;
  case 3:
    calc_g3_av(&gx, window, weights); break;
  default: ;
  }
  for (i=0; i<n_configs; i++) { 
    sprintf(buffer,"%f ",gx[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }
  free(gx);
  return (TCL_OK);
}


int tclcommand_analyze_parse_formfactor(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { formfactor | <formfactor> } <qmin> <qmax> <qbins> [<chain_start> <n_chains> <chain_length>]' */
  /***********************************************************************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+5];
  int i;
  double qmin,qmax, q,qfak, *ff; int qbins;
  if (argc < 3) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze formfactor <qmin> <qmax> <qbins> [<chain_start> <n_chains> <chain_length>]",
		     (char *)NULL);
    return (TCL_ERROR);
  } else {
    if (!ARG0_IS_D(qmin))
      return (TCL_ERROR);
    if (!ARG1_IS_D(qmax))
      return (TCL_ERROR);
    if (!ARG_IS_I(2, qbins))
      return (TCL_ERROR);
    argc-=3; argv+=3;
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;

  if ((chain_n_chains == 0) || (chain_length == 0)) {
    Tcl_AppendResult(interp, "The chain topology has not been set",(char *)NULL); return TCL_ERROR;
  }
  
  if (qbins <=0) {
    Tcl_AppendResult(interp, "Nothing to be done - choose <qbins> greater zero to get S(q)!",(char *)NULL); return TCL_ERROR;
  }

  if (qmin <= 0.) {
    Tcl_AppendResult(interp, "formfactor S(q) requires qmin > 0", (char *)NULL);
    return TCL_ERROR;
  }
  if (qmax <= qmin) {
    Tcl_AppendResult(interp, "formfactor S(q) requires qmin < qmax", (char *)NULL);
    return TCL_ERROR;
  }

  if (!average) analyze_formfactor(qmin, qmax, qbins, &ff);
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze formfactor ...' to only look at current state!",
		     (char *)NULL);
    return TCL_ERROR; }
  else analyze_formfactor_av(qmin, qmax, qbins, &ff);
  
  q = qmin; qfak = pow((qmax/qmin),(1.0/qbins));
  for(i=0; i<=qbins; i++) { sprintf(buffer,"{%f %f} ",q,ff[i]); q*=qfak; Tcl_AppendResult(interp, buffer, (char *)NULL); }
  free(ff);
  return (TCL_OK);
}

int tclcommand_analyze_parse_rdfchain(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze { rdfchain } <r_min> <r_max> <r_bins> [<chain_start> <n_chains> <chain_length>]' */
  /***********************************************************************************************************/
  char buffer[4*TCL_DOUBLE_SPACE+7];
  int i, r_bins;
  double r_min, r_max, *f1, *f2, *f3;
  double bin_width, r;
  if (argc < 3) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze rdfchain <r_min> <r_max> <r_bins> [<chain_start> <n_chains> <chain_length>]",
		     (char *)NULL);
    return (TCL_ERROR);
  } else {
    if (!ARG0_IS_D(r_min))
      return (TCL_ERROR);
    if (!ARG1_IS_D(r_max))
      return (TCL_ERROR);
    if (!ARG_IS_I(2, r_bins))
      return (TCL_ERROR);
    argc-=3; argv+=3;
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  
  if ((chain_n_chains == 0) || (chain_length == 0)) {
    Tcl_AppendResult(interp, "The chain topology has not been set",(char *)NULL); return TCL_ERROR;
  }
  
  if (r_bins <=0) {
    Tcl_AppendResult(interp, "Nothing to be done - choose <r_bins> greater zero!",(char *)NULL); return TCL_ERROR;
  }

  if (r_min <= 0.) {
    Tcl_AppendResult(interp, "<r_min> has to be positive", (char *)NULL);
    return TCL_ERROR;
  }
  if (r_max <= r_min) {
    Tcl_AppendResult(interp, "<r_max> has to be larger than <r_min>", (char *)NULL);
    return TCL_ERROR;
  }
  updatePartCfg(WITHOUT_BONDS);
  analyze_rdfchain(r_min, r_max, r_bins, &f1, &f2, &f3);

  bin_width = (r_max - r_min) / (double)r_bins;
  r = r_min + bin_width/2.0;
  for(i=0; i<r_bins; i++) {
    sprintf(buffer,"{%f %f %f %f} ",r,f1[i],f2[i],f3[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
     r+= bin_width;
  }  
  free(f1); free(f2); free(f3);
  return (TCL_OK);
}

#ifdef ELECTROSTATICS
int tclcommand_analyze_parse_cwvac(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze { cwvac } <maxtau> <interval> [<chain_start> <n_chains> <chain_length>]' */
  /***********************************************************************************************************/
  char buffer[4*TCL_DOUBLE_SPACE+7];
  int i, maxtau, interval;
  double *avac, *evac;
  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze gkmobility <maxtau> <interval> [<chain_start> <n_chains> <chain_length>]",
		     (char *)NULL);
    return (TCL_ERROR);
  } else {
    if (!ARG0_IS_I(maxtau))
      return (TCL_ERROR);
    if (!ARG1_IS_I(interval))
      return (TCL_ERROR);
    argc-=2; argv+=2;
  }
  if (tclcommand_analyze_set_parse_chain_topology_check(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  
  if ((chain_n_chains == 0) || (chain_length == 0)) {
    Tcl_AppendResult(interp, "The chain topology has not been set",(char *)NULL); return TCL_ERROR;
  }
  
  if (maxtau <=0) {
    Tcl_AppendResult(interp, "Nothing to be done - choose <maxtau> greater zero!",(char *)NULL); return TCL_ERROR;
  }

  if (interval <= 0) {
    Tcl_AppendResult(interp, "<interval> has to be positive", (char *)NULL);
    return TCL_ERROR;
  }
  updatePartCfg(WITHOUT_BONDS);
  analyze_cwvac(maxtau, interval, &avac, &evac);
  // create return string
  sprintf(buffer, "{ ");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  for(i=0;i<=maxtau;i++) {
    sprintf(buffer,"%e ",avac[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer, "} { ");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  for(i=0;i<=maxtau;i++) {
    sprintf(buffer,"%e ",evac[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer, "}");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  free(avac); free(evac);
  return (TCL_OK);
}
#endif

