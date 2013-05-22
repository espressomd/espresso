/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#include "utils.h"
#include "parser.h"
#include <strings.h>

static void setup_linear_bins(DoubleList *dl, double min_bin, double max_bin, int bins)
{
  int i;
  realloc_doublelist(dl, dl->n = bins + 1);
  for (i = 0; i <= bins; i++)
    dl->e[i] = min_bin + ((max_bin - min_bin)/bins)*i;
}

static void setup_log_bins(DoubleList *dl, double min_bin, double max_bin, int bins)
{
  int i;
  realloc_doublelist(dl, dl->n = bins + 1);
  for (i = 0; i <= bins; i++)
    dl->e[i] = min_bin*pow(max_bin/min_bin, ((double)i)/bins);
}

int tclcommand_bin(ClientData cdata, Tcl_Interp *interp,
	int argc, char **argv)
{
  DoubleList coords, data, count, sum, bins;
  int i, num_bins, give_bincounts = 0;
  double min_bin, max_bin, contr;
  int w, s, e, c;
  char buffer[2 + TCL_DOUBLE_SPACE];

  init_doublelist(&coords);
  init_doublelist(&data);
  init_doublelist(&sum);
  init_doublelist(&bins);

  /* type of the binning */
  if (argc > 1 && ARG1_IS_S("-stats")) {
    give_bincounts = 1;
    argc -= 1; argv += 1;
  }

  /* type of the binning */
  if (argc > 2 && ARG1_IS_S("-bins")) {
    if (!ARG_IS_DOUBLELIST(2, bins))
      return TCL_ERROR;
    argc -= 2; argv += 2;
  }
  else if (argc > 4 &&
	   (ARG1_IS_S("-linbins") || ARG1_IS_S("-logbins"))) {
    if (!ARG_IS_D(2, min_bin) ||
	!ARG_IS_D(3, max_bin) ||
	!ARG_IS_I(4, num_bins))
      return TCL_ERROR;    

    /* swap if necessary */
    if (min_bin > max_bin) { double tmp = min_bin; min_bin = max_bin; max_bin = tmp; }

    if (ARG1_IS_S("-linbins")) setup_linear_bins(&bins, min_bin, max_bin, num_bins);
    else if (ARG1_IS_S("-logbins")) setup_log_bins(&bins, min_bin, max_bin, num_bins);

    argc -= 4; argv += 4;
  }

  if (bins.n < 2) {
    Tcl_AppendResult(interp, "please specify at least two bin boundaries", (char *) NULL);
    return TCL_ERROR;
  }

  /* determine job to do */
  if (argc == 2 && ARG1_IS_S("-binctrwdth")) {
    /* just return the centers */
    Tcl_PrintDouble(interp, .5*(bins.e[0] + bins.e[1]), buffer);
    Tcl_AppendResult(interp, "{", buffer, (char *) NULL);
    Tcl_PrintDouble(interp, bins.e[1] - bins.e[0], buffer);
    Tcl_AppendResult(interp, " ", buffer, "}", (char *) NULL);
    for (i = 1; i < bins.n - 1; i++) {
      Tcl_PrintDouble(interp, .5*(bins.e[i] + bins.e[i+1]), buffer);
      Tcl_AppendResult(interp, " {", buffer, (char *) NULL);
      Tcl_PrintDouble(interp, bins.e[i+1] - bins.e[i], buffer);
      Tcl_AppendResult(interp, " ", buffer, "}", (char *) NULL);
    }
    return TCL_OK;
  }
  else if (argc == 2) {
    /* do the binning with bisection search algorithm */

    /* check for the type of data */
    if (!ARG1_IS_DOUBLELIST(coords)) {
      int i, tmp_argc, parse_error = 0;
      char  **tmp_argv;
      Tcl_ResetResult(interp);
      Tcl_SplitList(interp, argv[1], &tmp_argc, &tmp_argv);
      realloc_doublelist(&coords, coords.n = tmp_argc);
      realloc_doublelist(&data, data.n = tmp_argc);
      for(i = 0 ; i < tmp_argc; i++) {
	int tmp_argc2;
	char  **tmp_argv2;
	Tcl_SplitList(interp, tmp_argv[i], &tmp_argc2, &tmp_argv2);
	if (tmp_argc2 != 2) {
	  Tcl_AppendResult(interp, "data set has to be either a list of doubles or of lists of 2 doubles", (char *) NULL);
	  parse_error = 1; break;
	}
	if (Tcl_GetDouble(interp, tmp_argv2[0], &(coords.e[i])) == TCL_ERROR) { parse_error = 1; break; }
	if (Tcl_GetDouble(interp, tmp_argv2[1], &(data.e[i])) == TCL_ERROR) { parse_error = 1; break; }
	Tcl_Free((char *)tmp_argv2);
      }
      Tcl_Free((char *)tmp_argv);
      if (parse_error) return TCL_ERROR;
    }
      
    /* the binning itself */
    alloc_doublelist(&count, count.n = bins.n - 1);
    for (i = 0; i < count.n; i++) count.e[i] = 0;
    if (data.n) {
      alloc_doublelist(&sum, sum.n = bins.n - 1);
      for (i = 0; i < sum.n; i++) sum.e[i] = 0;
    }

    for (i = 0; i < coords.n; i++) {
      double cd = coords.e[i];
      if (cd < bins.e[0] || cd > bins.e[bins.n-1])
	continue;
      s = 0; e = bins.n - 1;
      while ((w = e - s) > 1) {
	c = (e + s)/2;
	if (cd >= bins.e[c]) s = c; else e = c;
      }
      count.e[s]++;
      if (data.n)
	sum.e[s] += data.e[i];
    }
    
    /* normalization */
    contr = 1./coords.n;

    for (i = 0; i < count.n; i++) {
      if (data.n) {
	if (count.e[i]) {
	  double tmp = sum.e[i]/count.e[i];
	  Tcl_PrintDouble(interp, tmp, buffer);
	}
	else
	  strcpy(buffer, "n/a");
      }
      else {
	Tcl_PrintDouble(interp, count.e[i] * contr, buffer);
      }
 
      if (i == 0)
	Tcl_AppendResult(interp, buffer, (char *) NULL);
      else
	Tcl_AppendResult(interp, " ", buffer, (char *) NULL);

      if (give_bincounts) {
	sprintf(buffer, "%d", (int)count.e[i]);
	Tcl_AppendResult(interp, " ", buffer, (char *) NULL);
      }
    }
    return TCL_OK;
  }

  Tcl_ResetResult(interp);
  Tcl_AppendResult(interp, "usage: bin -bins <binboundarylist> | "
		   "(-linbins|-logbins <start> <end> <num>) <data>|-binctrwdth\n", (char *) NULL);
  Tcl_AppendResult(interp, "       <data> is a list of doubles to bin or lists {coord data},"
		   " where data is to be averaged in each bin", (char *) NULL);
  return TCL_ERROR;   
}

