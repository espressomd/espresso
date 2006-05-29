// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file uwerr.c
    Implementation of \ref uwerr.h "uwerr.h".
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include "utils.h"
#include "uwerr.h"


/** This is eps from matlab */
#define UW_EPS 2.2204e-16

/*
enum UWerr_err_t {
  UW_NO_ERROR              = 0x000,
  UW_INVALID_DATA          = 0x001,
  UW_INVALID_ROWS          = 0x002,
  UW_INVALID_COLS          = 0x004,
  UW_INVALID_COL           = 0x008,
  UW_INVALID_FUNC          = 0x010,
  UW_INVALID_NREP          = 0x020,
  UW_INVALID_LEN           = 0x040,
  UW_NO_MEMORY             = 0x080,
  UW_WINDOWING_COND_FAILED = 0x100,
  UW_SIGMA_BIAS_CANCELLED  = 0x200
};
*/

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

/*
char * UWerr_err_to_str(struct UWerr_t s) {
  char * err_str[] = {"No error.",
		      "Invalid data.",
		      "Wrong number of rows for data.",
		      "Wrong number of columns for data.",
		      "Invalid function given.",
		      "Wrong vector for representations.",
		      "Wrong number for length of representations vector.",
		      "Out of memory.",
		      "Windowing condition failed at %d.",
		      "Sigma bias cancelled at %.3f."};

  return 0;
}
*/

double gammaln(double a)
{
    int j;
    double x, y, tmp, ser;
    double cof[6] = { 76.18009172947146,
                     -86.50532032941677,
                      24.01409824083091,
                      -1.231739572450155,
                       0.1208650973866179e-2,
                      -0.5395239384953e-5};

    y = x = a;
    tmp = x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j < 6; ++j)
	ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
}

void gammaser(double * gser, double a, double x)
{
    int n, ITMAX = 100;
    double sum, del, ap, gln, eps = DBL_EPSILON;

    gln = gammaln(a);

    if (x <= 0.0) {
	if (x < 0.0)
	    puts("uwerr: x less than 0 in gammaser.");
	*gser = 0.0;
	return;
    } else {
	ap = a;
	del = sum = 1.0/a;
	for (n = 0; n < ITMAX; ++n) {
	    ++ap;
	    del *= x/ap;
	    sum += del;
	    if (fabs(del) < fabs(sum)*eps) {
		*gser = sum*exp(-x+a*log(x)-gln);
		return;
	    }
	}
	puts("uwerr: a too large, ITMAX too small in gammaser.");
    }
}


void gammacf(double * gcf, double a, double x)
{
    int i, ITMAX = 100;
    double an, b, c, d, del, h, gln, eps = DBL_EPSILON, FPMIN = DBL_MIN;

    gln=gammaln(a);
    b = x+1.0-a;
    c = 1.0/FPMIN;
    d = 1.0/b;
    h = d;
    for (i = 1; i <= ITMAX; ++i) {
	an = -i*(i-a);
	b += 2.0;
	d = an*d+b;
	if (fabs(d) < FPMIN)
	    d = FPMIN;
	c = b+an/c;
	if (fabs(c) < FPMIN)
	    c = FPMIN;
	d = 1.0/d;
	del = d*c;
	h *= del;
	if (fabs(del-1.0) <= eps)
	    break;
    }
    if (i > ITMAX)
	puts("uwerr: a too large, ITMAX too small in gammacf.");
    *gcf = exp(-x+a*log(x)-gln)*h;
}

/** The incomplete Gammafunction.

    This is the implementation of the incomplete Gammafunction as described in
    Nummerical Recepies in C++ The Art of Scientific Computing Second Edition.

    The incomplete Gammafunction is defined as
    \f[
       Q(a,x):=\frac{1}{\Gamma(a)}\int_x^\infty e^{-t}t^{a-1}\textrm{d}t\quad(a>0)
    \f]
 */

double gammaq(double a, double x)
{
    double retval=0;

    if (x < 0.0 || a <= 0.0) {
	puts("uwerr: Invalid arguments for gammaq.");
	return 0.0;
    }

    if (x < a+1.0) {
	gammaser(&retval, a, x);
	return 1.0-retval;
    }
    /* else */
    gammacf(&retval, a, x);
    return retval;
}

/** Sum up a vector.

    \param v A pointer to the vector.
    \param len The length of the vector.
    \return The sum of all vector elements.
 */
double UWerr_sum(double * v, int len)
{
  int i;
  double s = 0;
  
  if (v)
    for (i = 0; i < len; ++i)
      s += v[i];

  return s;
}

/** Sum up the product of two vectors.

    \f$$\sum_{k=1}^\textrm{len}v[k]w[k]$\f$

    \param v A pointer to the first vector.
    \param w A pointer to the second vector.
    \param len The length of the shortest vector.
    \return The sum of the products of vector elements.
 */
double UWerr_dsum_double(double * v, double * w, int len)
{
  int i;
  double s = 0;
  
  if (v && w)
    for (i = 0; i < len; ++i)
      s += v[i]*w[i];

  return s;
}

/** Sum up the product of two vectors.

    \f$$\sum_{k=1}^\textrm{len}v[k]w[k]$\f$

    \param v A pointer to the first vector.
    \param w A pointer to the second vector.
    \param len The length of the shortest vector.
    \return The sum of the products of vector elements.
 */
double UWerr_dsum_int(int * v, double * w, int len)
{
  int i;
  double s = 0;
  
  if (v && w)
    for (i = 0; i < len; ++i)
      s += v[i]*w[i];

  return s;
}

/* forward declaration implementation below */
int uwerr_read_double_vector(Tcl_Interp *interp, char * data_in ,
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
int UWerr_proj(ClientData cd, Tcl_Interp *interp, int argc, char *argv[])
{
  double * a;
  int len, n;
  char ret[TCL_DOUBLE_SPACE];

  if (argc != 3)
    return TCL_ERROR;

  if (uwerr_read_double_vector(interp, argv[1], &a, &len) == TCL_ERROR) {
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
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, ret, (char *)NULL);
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
int UWerr_f(Tcl_Interp *interp, Tcl_CmdInfo * cmdInfo, int argc, char ** argv,
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
  char ** my_argv;

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

  if (!(my_argv=(char**)malloc((argc+1)*sizeof(char*)))) {
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

  Fbb = strtod(interp->result,0);
  for (k = 0; k < len; ++k) {
    uwerr_write_tcl_vector(interp, abr[k], cols, tcl_vector);
    Tcl_ResetResult(interp);

    if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
      goto err_exit;
    Fbr[k] = strtod(interp->result,0);
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
      fgrad[a] = strtod(interp->result,0);

      abb[a] = tmp - std_a;

      uwerr_write_tcl_vector(interp, abb, cols, tcl_vector);
      Tcl_ResetResult(interp);
      if (cmdInfo->proc(cmdInfo->clientData, interp, argc+1, my_argv) != TCL_OK)
	goto err_exit;
      fgrad[a] -= strtod(interp->result,0);

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
  char * argv[2];
  char * name = "UWerrInternalFunction";
  int res;
  
  argv[0] = name;
  argv[1] = (char*)malloc(TCL_INTEGER_SPACE*sizeof(char));
  sprintf(argv[1], "%d", col_to_analyze);

  Tcl_CreateCommand(interp, name, UWerr_proj, 0, NULL);
  Tcl_GetCommandInfo(interp, name, &cmdInfo);
  
  res = UWerr_f(interp, &cmdInfo, 2, argv,
		data, rows, cols, n_rep, len, s_tau, plot);

  Tcl_DeleteCommand(interp, name);
  
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

  if (Tcl_SplitList(interp, data_in, nrows, &row) == TCL_ERROR)
    return TCL_ERROR;

  if (*nrows < 1) {
    Tcl_AppendResult(interp, "first argument has to be a matrix.",
		     (char *)NULL);
    return TCL_ERROR;
  }

  if (!(*data = (double**)malloc(*nrows*sizeof(double*)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)NULL);
    Tcl_Free((char *)row);
    return TCL_ERROR;
  }

  for (i = 0; i < *nrows; ++i) {
    tmp_ncols = -1;
    
    if (Tcl_SplitList(interp, row[i], &tmp_ncols, &col) == TCL_ERROR) {
      Tcl_Free((char*)row);
      return TCL_ERROR;
    }

    if (i == 0) {
      if (tmp_ncols < 1) {
	Tcl_AppendResult(interp, "first argument has to be a matrix.",
			 (char *)NULL);
	Tcl_Free((char *)col);
	Tcl_Free((char*)row);
	return TCL_ERROR;
      }

      *ncols = tmp_ncols;

    } else if (*ncols != tmp_ncols) {
      Tcl_AppendResult(interp, "number of columns changed.",
		       (char *)NULL);
      Tcl_Free((char *)col);
      Tcl_Free((char*)row);
      return TCL_ERROR;
    }

    if (!((*data)[i] = (double*)malloc(*ncols*sizeof(double)))) {
      Tcl_AppendResult(interp,"Out of Memory.",
		       (char *)NULL);
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

  if (Tcl_SplitList(interp, data_in, len, &col) == TCL_ERROR)
    return TCL_ERROR;

  if (*len < 1) {
    Tcl_AppendResult(interp, "Argument is not a vector.",
		     (char *)NULL);
    return TCL_ERROR;
  }

  if (!(*nrep = (int*)malloc((*len)*sizeof(int)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)NULL);
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
int uwerr_read_double_vector(Tcl_Interp *interp, char * data_in ,
			     double ** nrep, int * len)
{
  char ** col;
  int i;

  *len = -1;

  if (Tcl_SplitList(interp, data_in, len, &col) == TCL_ERROR)
    return TCL_ERROR;

  if (*len < 1) {
    Tcl_AppendResult(interp, "Argument is not a vector.",
		     (char *)NULL);
    return TCL_ERROR;
  }

  if (!(*nrep = (double*)malloc((*len)*sizeof(double)))) {
    Tcl_AppendResult(interp, "Out of Memory.",
		     (char *)NULL);
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

int uwerr(ClientData cd, Tcl_Interp *interp, int argc, char *argv[])
{
  int i, nrows, ncols, len, plot = 0,
    col_to_analyze = -1, analyze_col = 0, error = 0,
    result = TCL_OK;
  double s_tau = 1.5;
  int * nrep;
  double ** data;
  char * str;
  char ** my_argv;
  Tcl_CmdInfo cmdInfo;

  if (argc < 4) {
    Tcl_AppendResult(interp, argv[0], " needs at least 3 arguments.\n",
		     "usage: ", argv[0], " <data> <nrep> {<col>|<f>} [<s_tau> [<f_args>]] [plot]\n",
		     (char *)NULL);
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
      error = 1;
      str = (char *)malloc(TCL_INTEGER_SPACE*sizeof(char));
      sprintf(str, "%d", ncols);
      Tcl_AppendResult(interp, "third argument has to be a function or a ",
		       "number between 1 and ", str, "!", (char *)NULL);
      free(str);
    }
  }
  
  if (!error && analyze_col &&
      (col_to_analyze < 1 || col_to_analyze > ncols)) {
    error = 1;
    str = (char *)malloc(TCL_INTEGER_SPACE*sizeof(char));
    sprintf(str, "%d", ncols);
    Tcl_AppendResult(interp, "third argument has to be a function or a ",
		     "number between 1 and ", str, ".", (char *)NULL);
    free(str);
  }

  /* check for plot as fourth argument */
  if (argc > 4 && !error) {
    if (!strcmp(argv[4], "plot"))
      plot = 1;
    else {

      /* read s_tau if there is a fourth arg */
      if (Tcl_GetDouble(interp, argv[4], &s_tau) == TCL_ERROR) {
	error = 1;
	Tcl_AppendResult(interp, "fourth argument has to be a double or 'plot'.", (char *)NULL);
      }

    }
  }

  if (argc > 5 && ! error)
    if (!strcmp(argv[argc-1], "plot"))
      plot = 1;

  if (!error && analyze_col) {
    result = UWerr(interp, data, nrows, ncols,
		   col_to_analyze-1, nrep, len, s_tau, plot);
  }

  if (!error && !analyze_col) {
    my_argv = (char**)malloc((argc-3)*sizeof(char*));
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

  return error ? TCL_ERROR : TCL_OK;
}
