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
/** \file uwerr_tcl.cpp
    Implements the uwerr command.
*/

#include <cstdlib>
#include <cstring>
#include <cmath>
#include "uwerr.hpp"
#include "config.hpp"
#include <tcl.h>


/** This is eps from matlab */
#define UW_EPS 2.2204e-16

struct UWerr_t {
  double value;
  double dvalue;
  double ddvalue;
  double tau_int;
  double dtau_int;
  double Q_val;
  int error;
  int W;
  double bias;
};

/* forward declaration implementation below */
int uwerr_read_tcl_double_vector(Tcl_Interp *interp, const char * data_in ,
			     double ** data_out, int * len);

/** Create a string with enough space for in_len doubles
 */
int uwerr_create_tcl_vector(char ** out, int in_len)
{
  *out = (char*)malloc(in_len*(TCL_DOUBLE_SPACE+1)*sizeof(char)+1);
  if (!(*out))
    return -1;

  return 0;
}

/** Write a vector of doubles as string
 */
int uwerr_write_tcl_vector(Tcl_Interp * interp,
			   double * in, int in_len, char * out)
{
  char str[TCL_DOUBLE_SPACE];
  int i, cur_len = 0;

  if (!out)
    return -1;

  out[0] = 0;

  for (i = 0; i < in_len; ++i) {
    Tcl_PrintDouble(interp, in[i], str);
    strcat(out+cur_len, str);
    cur_len += strlen(str);
    out[cur_len] = ' '; // overwrite the trailing NULL
    out[cur_len+1] = 0;
    cur_len += 1;
  }

  return 0;
}

/** Free the resources for the given string
 */
int uwerr_free_tcl_vector(char * vec) {
  if (!vec)
    return -1;
  
  free(vec);
  return 0;
}

/** Projection of a vector to one entry.

    This function is a replacement for a tcl function like
    \verbatim proc proj {vec n} { return [lindex $vec $n] }\endverbatim

    Therefor \em argc has to be 3 and \em argv has to consist of
    <ul><li>a name of the function (e.g. in \ref UWerr)
        <li>the vector as Tcl list
        <li>the coordinate to project to</ul>
 */
int UWerr_proj(ClientData cd, Tcl_Interp *interp, int argc, const char *argv[])
{
  double * a;
  int len, n;
  char ret[TCL_DOUBLE_SPACE];

  if (argc != 3)
    return TCL_ERROR;

  if (uwerr_read_tcl_double_vector(interp, argv[1], &a, &len) == TCL_ERROR) {
    Tcl_AppendResult(interp, "\nFirst arg to UWerr_proj is wrong.", (char *)NULL);
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &n) == TCL_ERROR) {
    free(a);
    Tcl_AppendResult(interp, "Second arg to UWerr_proj is wrong.", (char *)NULL);
    return TCL_ERROR;
  }

  if (n < len && n >= 0) {
    Tcl_PrintDouble(interp, a[n], ret);
    Tcl_SetResult(interp, ret, TCL_VOLATILE);
    return TCL_OK;
  }

  free(a);
  Tcl_AppendResult(interp, "Second arg to UWerr_proj is out of range", (char *)NULL);
  return TCL_ERROR;
}

/** The main function.

    The function implementing the algorithm described in

    arXiv:hep-lat/0306017 v1 13 Jun 2003 \em Wolff, U. \em Monte Carlo errors with less errors.
*/
int UWerr_f(Tcl_Interp *interp, Tcl_CmdInfo * cmdInfo, int argc, const char ** argv,
	    double ** data, int rows, int cols,
	    int * n_rep, int len, double s_tau, int plot)
{
  struct UWerr_t ret;
  int a, k, i, sum = 0, W_opt = 0, W_max = 0;
  double Fbb = 0, bF = 0, Fb = 0, * abb = 0L, tau = 0, tmp;
  double ** abr = 0L, * Fbr = 0L, * fgrad = 0L, * delpro = 0L;
  double * gFbb = 0L, CFbb_opt = 0, G_int = 0, std_a;
  char flag = 0;
  char * str = 0L;
  char * tcl_vector = 0L;
  const char ** my_argv = (const char**)malloc((argc+1)*sizeof(char*));

  FILE * plotDataf, * plotScriptf;

  ret.Q_val = 0;

  if (!data) {
    Tcl_AppendElement(interp, "No data matrix given.");
    return TCL_ERROR;
  }
  if (rows < 1) {
    Tcl_AppendElement(interp, "Data matrix has no rows.");
    return TCL_ERROR;
  }
  if (cols < 1) {
    Tcl_AppendElement(interp, "Data matrix has no columns.");
    return TCL_ERROR;
  }
  if(!cmdInfo && !cmdInfo->proc) {
    Tcl_AppendElement(interp, "No function to call given.");
    return TCL_ERROR;
  }
  if (!n_rep) {
    Tcl_AppendElement(interp, "No representations vector given.");
    return TCL_ERROR;
  }
  if (len < 1) {
    Tcl_AppendElement(interp, "Representations vector is empty.");
    return TCL_ERROR;
  }
  
  /* \sum_{i=1}^{len} n_rep[i-1] = rows */

  k = rows; /* for now k is going to be min(n_rep) */
  for (i = 0; i < len; ++i) {
    sum += n_rep[i];
    if (n_rep[i] < k)
      k = n_rep[i];
  }

  if (sum != rows || k <= 0) {
    Tcl_AppendElement(interp, "Representations vector is invalid.");
    return TCL_ERROR;
  }

  if (s_tau > 0) {
    W_max = (int)rint(k/2.); /* until here: k = min(n_rep) */
    flag = 1;
    if (W_max < 1) W_max = 1;
  }

  /* string for output of numbers */
  str = (char *)malloc((TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE)*sizeof(char));

  if (!(delpro = (double*)malloc(rows*sizeof(double)))) {
    Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
    free(str);
    return TCL_ERROR;
  }

  if (!(Fbr = (double*)malloc(len*sizeof(double)))) {
    free(delpro);
    free(str);
    Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
    return TCL_ERROR;
  }

  if (!(fgrad = (double*)malloc(cols*sizeof(double)))) {
    free(delpro);
    free(Fbr);
    free(str);
    Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
    return TCL_ERROR;
  }

  if (!(abb = (double*)malloc(cols*sizeof(double)))) {
    free(delpro);
    free(Fbr);
    free(fgrad);
    free(str);
    Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
    return TCL_ERROR;
  }

  /* abr \in (\Real)_{len, cols} */
  if (!(abr = (double**)malloc(len*sizeof(double*)))) {
    free(delpro);
    free(Fbr);
    free(fgrad);
    free(abb);
    free(str);
    Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
    return TCL_ERROR;
  }
  for (i = 0; i < len; ++i)
    if (!(abr[i] = (double*)malloc(cols*sizeof(double)))) {
      for (k = 0; k < i; ++k)
	free(abr[k]);

      free(abr);
      free(delpro);
      free(Fbr);
      free(fgrad);
      free(abb);
      free(str);
    
      Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
      return TCL_ERROR;
    }

  
  if (W_max > 0) {
    if (!(gFbb = (double*)malloc((W_max+1)*sizeof(double)))) {
      free(delpro);
      free(Fbr);
      free(fgrad);
      free(abb);
      for (k = 0; k < len; ++k)
	free(abr[k]);
      free(abr);

      free(str);
      Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
      return TCL_ERROR;
    }
  }
  
  if (uwerr_create_tcl_vector(&tcl_vector, cols)) {
      free(delpro);
      free(Fbr);
      free(fgrad);
      free(abb);
      for (k = 0; k < len; ++k)
	free(abr[k]);
      free(abr);
      free(gFbb);

      free(str);
      Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
      return TCL_ERROR;
  }

  if (!my_argv) {
      free(delpro);
      free(Fbr);
      free(fgrad);
      free(abb);
      for (k = 0; k < len; ++k)
	free(abr[k]);
      free(abr);
      free(gFbb);

      free(str);
      uwerr_free_tcl_vector(tcl_vector);
      Tcl_AppendResult(interp, "Out of Memory.", __LINE__, (char *) NULL);
      return TCL_ERROR;
  }

  my_argv[0] = argv[0];
  my_argv[1] = tcl_vector;
  for (i = 1; i < argc; ++i)
    my_argv[i+1] = argv[i];


  /* first we calculate N_r\bar{a}_\alpha^r \forall r, alpha */
  
  sum = 0;
  for (k = 0; k < len; ++k) {
    for (i = 0; i < n_rep[k]; ++i) {
      for (a = 0; a < cols; ++a) {
	if (i > 0)
	  abr[k][a] += data[sum + i][a];
	else
	  abr[k][a] = data[sum][a];
      }
    }
    sum += n_rep[k];
  }

  /* now we calculate \bar{\bar{a}}_\alpha \forall \alpha */

  for (k = 0; k < len; ++k) {
    for (a = 0; a < cols; ++a) {
      if (k > 0)
	abb[a] += abr[k][a];
      else
	abb[a] = abr[k][a];
    }
  }
  for (a =0; a < cols; ++a)
    abb[a] /= rows;

  /* now we calculate \bar{a}_\alpha^r with \forall \alpha */
  for (k = 0; k < len; ++k)
    for (a = 0; a < cols; ++a)
      abr[k][a] /= n_rep[k];

  uwerr_write_tcl_vector(interp, abb, cols, tcl_vector);
  Tcl_ResetResult(interp);

  if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
    goto err_exit;

  Fbb = strtod(Tcl_GetStringResult(interp),0);
  for (k = 0; k < len; ++k) {
    uwerr_write_tcl_vector(interp, abr[k], cols, tcl_vector);
    Tcl_ResetResult(interp);

    if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
      goto err_exit;
    Fbr[k] = strtod(Tcl_GetStringResult(interp),0);
  }

  Fb  = UWerr_dsum_int(n_rep, Fbr, len);
  Fb /= rows;

  for (a = 0; a < cols; ++a) {
    std_a = 0;
    for (k = 0; k < rows; ++k)
      std_a += (data[k][a]-abb[a])*(data[k][a]-abb[a]);
    std_a = sqrt(std_a)/rows;

    
    /* calc the gradient of f using df/da ~ (f(a+h)-f(a-h))/2*h 
       where h is the standard deviation divided by the sqrt of the 
       number of samples (= rows).
       Remember: abb[a] is the average for column a of data */

    if (std_a == 0)
      fgrad[a] = 0;
    else {
      tmp = abb[a];
      abb[a] += std_a;

      uwerr_write_tcl_vector(interp, abb, cols, tcl_vector);
      Tcl_ResetResult(interp);
      if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
	goto err_exit;
      fgrad[a] = strtod(Tcl_GetStringResult(interp),0);

      abb[a] = tmp - std_a;

      uwerr_write_tcl_vector(interp, abb, cols, tcl_vector);
      Tcl_ResetResult(interp);
      if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
	goto err_exit;
      fgrad[a] -= strtod(Tcl_GetStringResult(interp),0);

      abb[a] = tmp;
      fgrad[a] /= 2*std_a;
    }
  }

  /* calc delpro = data*fgrad - abb.*fgrad and
     the mean of delpro.^2 = gFbb[0] */

  tmp = UWerr_dsum_double(abb, fgrad, cols);
  gFbb[0] = 0;
  for (i = 0; i < rows; ++i) {
    delpro[i] = 0;

    for (a = 0; a < cols; a++) {
      delpro[i] += data[i][a]*fgrad[a];
    }
    delpro[i] -= tmp;

    gFbb[0] += delpro[i]*delpro[i];
  }
  gFbb[0] /= rows;

  i = 0;
  while(i < W_max) {
    gFbb[i+1] = 0;
    sum = 0;
    for (k = 0; k < len; ++k) {
      gFbb[i+1] += UWerr_dsum_double(delpro + sum, delpro + sum + i + 1, n_rep[k]-i-1);
      sum += n_rep[k];
    }
    gFbb[i+1] /= rows-(i+1)*len;

    if (flag) {
      G_int += gFbb[i+1]/gFbb[0];
      if (G_int <= 0)
	tau = UW_EPS;
      else
	tau = s_tau/log((G_int+1)/G_int);
      if (exp(-(i+1)/tau)-tau/sqrt((i+1)*rows) < 0) {
	W_opt = i+1;
	W_max = (W_max < 2*W_opt) ? W_max : 2*W_opt;
	flag = 0;
      }
    }
    ++i;
  }
  --i;

  if (flag) {
    W_opt = W_max;
    sprintf(str, "%d", W_max);
    Tcl_AppendResult(interp, "Windowing condition failed up to W = ", str, ".\n", (char *)NULL);
  }
  ret.W = W_opt;

  CFbb_opt = (gFbb[0] + 2*UWerr_sum(gFbb+1, W_opt))/rows;
  for (k = 0; k < i; ++k)
    gFbb[k] += CFbb_opt;
  CFbb_opt = (gFbb[0] + 2*UWerr_sum(gFbb+1, W_opt));

  ret.dvalue = sqrt(CFbb_opt/rows); /* sigmaF */
  
  if (len >= 2) {
    bF = (Fb-Fbb)/(len-1);
    Fbb -= bF;
    if (fabs(bF) > ret.dvalue/4) {
      Tcl_PrintDouble(interp, bF/ret.dvalue, str);
      Tcl_AppendResult(interp, "A ", str, " sigma bias of the mean has been cancelled./n", (char *)NULL);
    }
    for (i = 0; i < len; ++i)
      Fbr[i] -= bF*rows/n_rep[i];
    Fb -= bF*len;
    
    ret.bias = bF/ret.dvalue;
  }

  ret.tau_int = 0;
  for (i = 0; i <= W_opt; ++i)
    ret.tau_int += gFbb[i];
      
  ret.tau_int /= gFbb[0];
  ret.tau_int -= .5;

  ret.value  = Fbb;
  ret.ddvalue = ret.dvalue*sqrt((W_opt + .5)/rows);
  ret.dtau_int = 2 * ret.tau_int * sqrt((W_opt + .5 - ret.tau_int)/rows);

  if (len > 1) {
    for (i = 0; i < len; ++i)
      Fbr[i] = (Fbr[i] - Fb)*(Fbr[i] - Fb)*n_rep[i];
    
    ret.Q_val = UWerr_sum(Fbr, len);
    ret.Q_val /= CFbb_opt;
    ret.Q_val = gammaq((len-1)/2., ret.Q_val/2.);
  }

  if (plot) {
    plotScriptf = fopen("uwerr_plot_script", "w");

    fprintf(plotScriptf, "set ylabel \"Gamma\"; set xlabel \"W\"; set label \"W_opt=%d\" at %d,0 center; plot f(x) = 0, f(x) notitle, 'uwerr_plot_data' using 1:2 title \"normalized autocorrelation\" with lines; show label; pause -1\n", W_opt, W_opt);
    fprintf(plotScriptf, "set ylabel \"tau_int\"; plot f(x) = %.3f, 'uwerr_plot_data' using 1:3 title \"tau_int with statistical errors\" with lines,", ret.tau_int);
    fprintf(plotScriptf, " 'uwerr_plot_data' using 1:3:4 notitle with errorbars, f(x) title \"estimate\"; pause -1\n");

    fclose(plotScriptf);

    plotDataf = fopen("uwerr_plot_data", "w");
    tmp = 0;
    for (i = 0; i < W_max; ++i) {
      tmp += gFbb[i];
      /* print values for x-Axis, Gamma/Gamma[0], tau_int, and its errors */
      fprintf(plotDataf, "%d %.3f %.3f %.3f\n", i, gFbb[i]/gFbb[0],
	      tmp/gFbb[0]-.5, 2*sqrt((i+tmp/gFbb[0])/rows));
    }
    fclose(plotDataf);

    puts("Press Return to continue ...");
    Tcl_Eval(interp, "[exec gnuplot uwerr_plot_script]");
  }

  Tcl_ResetResult(interp);
  Tcl_PrintDouble(interp, ret.value, str);
  Tcl_AppendResult(interp, str, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ret.dvalue, str);
  Tcl_AppendResult(interp, str, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ret.ddvalue, str);
  Tcl_AppendResult(interp, str, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ret.tau_int, str);
  Tcl_AppendResult(interp, str, " ", (char *)NULL);
  Tcl_PrintDouble(interp, ret.dtau_int, str);
  Tcl_AppendResult(interp, str, (char *)NULL);
  if (len > 1) {
    Tcl_PrintDouble(interp, ret.Q_val, str);
    Tcl_AppendResult(interp, " ", str, (char *)NULL);
  }

 err_exit:
  free(abb);
  for (k = 0; k < len; ++k)
    free(abr[k]);
  free(abr);
  free(delpro);
  free(gFbb);
  free(Fbr);
  free(fgrad);
  free(str);
  free(my_argv);
  uwerr_free_tcl_vector(tcl_vector);

  return TCL_OK;
}

/** The function for analyzing a column only.
    \anchor UWerr
*/
int UWerr(Tcl_Interp * interp,
	  double ** data, int rows, int cols,
	  int col_to_analyze,
	  int * n_rep, int len,
	  double s_tau, int plot)
{
  Tcl_CmdInfo cmdInfo;
  char* argv[2];
  char* name = strdup("UWerrInternalFunction");
  int res;
  
  argv[0] = name;
  argv[1] = (char*)malloc(TCL_INTEGER_SPACE*sizeof(char));
  sprintf(argv[1], "%d", col_to_analyze);

  if (Tcl_CreateCommand(interp, name, UWerr_proj, 0, 0) == 0) {
      Tcl_AppendResult(interp, "could not create command \"", name, "\"", (char *)0);
      return TCL_ERROR;
  }
  if (Tcl_GetCommandInfo(interp, name, &cmdInfo) == 0) {
      Tcl_AppendResult(interp, "could not access command \"", name, "\"", (char *)0);
      return TCL_ERROR;
  }

  res = UWerr_f(interp, &cmdInfo, 2, (const char **)argv,
		data, rows, cols, n_rep, len, s_tau, plot);

  Tcl_DeleteCommand(interp, name);
  
  free(argv[0]);
  free(argv[1]);

  return res;
}

/** Reads a Tcl matrix and returns a C matrix.

    \param interp The Tcl interpreter
    \param data_in String containing a Tcl matrix of doubles
    \param data Pointer to the C matrix
    \param nrows Pointer to an int to store the height of the matrix
    \param ncols Pointer to an int to store the width of the matrix
    \return \em TCL_OK if everything went fine \em TCL_ERROR otherwise and 
            interp->result is set to an error message.

	    If \em TCL_OK is returned you have to make sure to free the memory
	    pointed to by data.
 */
int uwerr_read_matrix(Tcl_Interp *interp, char * data_in ,
		      double *** data, int * nrows, int * ncols)
{
  char ** row;
  char ** col;
  int tmp_ncols = -1, i, j, k;

  *nrows = *ncols = -1;

  if (Tcl_SplitList(interp, data_in, nrows, (const char ***)&row) == TCL_ERROR)
    return TCL_ERROR;

  if (*nrows < 1) {
    Tcl_AppendResult(interp, "first argument has to be a matrix.",
		     (char *)0);
    return TCL_ERROR;
  }

  if (!(*data = (double**)malloc(*nrows*sizeof(double*)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)0);
    Tcl_Free((char *)row);
    return TCL_ERROR;
  }

  for (i = 0; i < *nrows; ++i) {
    tmp_ncols = -1;
    
    if (Tcl_SplitList(interp, row[i], &tmp_ncols, (const char ***)&col) == TCL_ERROR) {
      Tcl_Free((char*)row);
      return TCL_ERROR;
    }

    if (i == 0) {
      if (tmp_ncols < 1) {
	Tcl_AppendResult(interp, "first argument has to be a matrix.",
			 (char *)0);
	Tcl_Free((char *)col);
	Tcl_Free((char*)row);
	return TCL_ERROR;
      }

      *ncols = tmp_ncols;

    } else if (*ncols != tmp_ncols) {
      Tcl_AppendResult(interp, "number of columns changed.",
		       (char *)0);
      Tcl_Free((char *)col);
      Tcl_Free((char*)row);
      return TCL_ERROR;
    }

    if (!((*data)[i] = (double*)malloc(*ncols*sizeof(double)))) {
      Tcl_AppendResult(interp,"Out of Memory.",
		       (char *)0);
      Tcl_Free((char *)row);
      Tcl_Free((char *)col);
      for (k = 0; k < i; ++k)
	free((*data)[i]);
      free(*data);
      return TCL_ERROR;
    };

    for (j = 0; j < *ncols; ++j) {
      if (Tcl_GetDouble(interp, col[j], &((*data)[i][j])) == TCL_ERROR) {
	Tcl_Free((char *)col);
	Tcl_Free((char *)row);
	for (k = 0; k <= i; ++k)
	  free((*data)[i]);
	free(*data);
	return TCL_ERROR;
      }
    }

    Tcl_Free((char *)col);
  }

  Tcl_Free((char *)row);

  return TCL_OK;
}

/** Reads a Tcl vector and returns a C vector.

    \param interp The Tcl interpreter
    \param data_in String containing a Tcl vector of integers
    \param nrep Pointer to the C vector
    \param len Pointer to an int to store the length of the vector
    \return \em TCL_OK if everything went fine \em TCL_ERROR otherwise and 
            interp->result is set to an error message.

	    If \em TCL_OK is returned you have to make sure to free the memory
	    pointed to by nrep.
 */
int uwerr_read_int_vector(Tcl_Interp *interp, char * data_in ,
			  int ** nrep, int * len)
{
  char ** col;
  int i;

  *len  = -1;
  *nrep =  0;

  if (Tcl_SplitList(interp, data_in, len, (const char ***)&col) == TCL_ERROR)
    return TCL_ERROR;

  if (*len < 1) {
    Tcl_AppendResult(interp, "Argument is not a vector.",
		     (char *)0);
    return TCL_ERROR;
  }

  if (!(*nrep = (int*)malloc((*len)*sizeof(int)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)0);
    Tcl_Free((char *)col);
    return TCL_ERROR;
  }

  for (i = 0; i < *len; ++i) {
      if (Tcl_GetInt(interp, col[i], &((*nrep)[i])) == TCL_ERROR) {
	Tcl_Free((char *)col);
	free(*nrep);
	return TCL_ERROR;
      }
  }

  Tcl_Free((char *)col);
  return TCL_OK;
}

/** Reads a Tcl vector and returns a C vector.

    \param interp The Tcl interpreter
    \param data_in String containing a Tcl vector of doubles
    \param nrep Pointer to the C vector
    \param len Pointer to an int to store the length of the vector
    \return \em TCL_OK if everything went fine \em TCL_ERROR otherwise and 
            interp->result is set to an error message.

	    If \em TCL_OK is returned you have to make sure to free the memory
	    pointed to by nrep.
 */
int uwerr_read_tcl_double_vector(Tcl_Interp *interp, const char * data_in ,
			     double ** nrep, int * len)
{
  char ** col;
  int i;

  *len = -1;

  if (Tcl_SplitList(interp, data_in, len, (const char ***)&col) == TCL_ERROR)
    return TCL_ERROR;

  if (*len < 1) {
    Tcl_AppendResult(interp, "Argument is not a vector.",
		     (char *)0);
    return TCL_ERROR;
  }

  if (!(*nrep = (double*)malloc((*len)*sizeof(double)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)0);
    Tcl_Free((char *)col);
    return TCL_ERROR;
  }

  for (i = 0; i < *len; ++i) {
      if (Tcl_GetDouble(interp, col[i], &((*nrep)[i])) == TCL_ERROR) {
	Tcl_Free((char *)col);
	free(*nrep);
	return TCL_ERROR;
      }
  }

  Tcl_Free((char *)col);
  return TCL_OK;
}

int tclcommand_uwerr(ClientData cd, Tcl_Interp *interp, int argc, char *argv[])
{
  int i, nrows, ncols, len, plot = 0,
    col_to_analyze = -1, analyze_col = 0,
    result = TCL_OK;
  double s_tau = 1.5;
  int * nrep;
  double ** data;
  char * str;
  Tcl_CmdInfo cmdInfo;

  if (argc < 4) {
    Tcl_AppendResult(interp, argv[0], " needs at least 3 arguments.\n",
		     "usage: ", argv[0], " <data> <nrep> {<col>|<f>} [<s_tau> [<f_args>]] [plot]\n",
		     (char *)0);
    return TCL_ERROR;
  }

  /* read the matrix containing the data */
  if (uwerr_read_matrix(interp, argv[1], &data, &nrows, &ncols) == TCL_ERROR)
    return TCL_ERROR;

  /* read the vector containing the length of each representation */
  if (uwerr_read_int_vector(interp, argv[2], &nrep, &len) == TCL_ERROR)
    return TCL_ERROR;

  /* check if we analyze a column or a function of the columns */
  if (!Tcl_GetCommandInfo(interp, argv[3], &cmdInfo)) {
    analyze_col = 1;
    if (Tcl_GetInt(interp, argv[3], &col_to_analyze) == TCL_ERROR) {
      result = TCL_ERROR;
      str = (char *)malloc(TCL_INTEGER_SPACE*sizeof(char));
      sprintf(str, "%d", ncols);
      Tcl_AppendResult(interp, "third argument has to be a function or a ",
		       "number between 1 and ", str, "!", (char *)0);
      free(str);
    }
  }
  
  if ((result == TCL_OK) && analyze_col &&
      (col_to_analyze < 1 || col_to_analyze > ncols)) {
    result = TCL_ERROR;
    str = (char *)malloc(TCL_INTEGER_SPACE*sizeof(char));
    sprintf(str, "%d", ncols);
    Tcl_AppendResult(interp, "third argument has to be a function or a ",
		     "number between 1 and ", str, ".", (char *)0);
    free(str);
  }

  /* check for plot as fourth argument */
  if (argc > 4 && (result != TCL_OK)) {
    if (!strcmp(argv[4], "plot"))
      plot = 1;
    else {

      /* read s_tau if there is a fourth arg */
      if (Tcl_GetDouble(interp, argv[4], &s_tau) == TCL_ERROR) {
	result = TCL_ERROR;
	Tcl_AppendResult(interp, "fourth argument has to be a double or 'plot'.", (char *)0);
      }

    }
  }

  if (argc > 5 && (result == TCL_OK))
    if (!strcmp(argv[argc-1], "plot"))
      plot = 1;

  if ((result == TCL_OK) && analyze_col) {
    result = UWerr(interp, data, nrows, ncols,
		   col_to_analyze-1, nrep, len, s_tau, plot);
  }

  if ((result == TCL_OK) && !analyze_col) {
      const char ** my_argv = (const char**)malloc((argc-3)*sizeof(char*));
    my_argv[0] = argv[3];
    for (i = 0; i < argc-5-plot; ++i)
      my_argv[i+1] = argv[5+i];
    result = UWerr_f(interp, &cmdInfo, argc-plot>4?argc-plot-4:1, my_argv,
		     data, nrows, ncols, nrep, len, s_tau, plot);
    free(my_argv);
  }
  
  for (i = 0; i < nrows; ++i)
    free(data[i]);
  free(data);

  free(nrep);

  return result;
}
