// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file specfunc.h
    This file contains implementations for some special functions which are needed by the MMM family of
    algorithms. This are the modified Hurwitz zeta function and the modified Bessel functions of first
    and second kind. The implementations are based on the GSL code (see \ref specfunc.c "specfunc.c"
    for the original GSL header).

    The Hurwitz zeta function is evaluated using the Euler-MacLaurin summation formula, the Bessel functions
    are evaluated using several different Chebychev expansions. Both achieve a precision of nearly machine
    precision, which is no problem for the Hurwitz zeta function, which is only used when determining the
    coefficients for the modified polygamma functions (see \ref mmm-common.h "mmm-common.h"). However, the
    Bessel functions are actually used in the near formula of MMM2D, which is therefore slightly slower than
    necessary. On the other hand, the number of terms in the Bessel sum is quite small normally, so that a less
    precise version will probably not generate a huge computational speed improvement.
*/
#ifndef SPECFUNC_H
#define SPECFUNC_H

/** Hurwitz zeta function. This function was taken from the GSL code. */
double hzeta(double order, double x);

/** Modified Bessel function of first kind, order 0. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
double I0(double x);
/** Modified Bessel function of first kind, order 1. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
double I1(double x);
/** Modified Bessel function of second kind, order 0. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
double K0(double x);
/** Modified Bessel function of second kind, order 1. This function was taken from
    the GSL code. Precise roughly up to machine precision. */
double K1(double x);

/** Besselfunctions K0 at x.
    The implementation has an absolute precision of around 10^(-14), which is
    comparable to the relative precision sqrt implementation of current hardware.
*/
double LPK0(double x);

/** Besselfunctions K1 at x.
    The implementation has an absolute precision of around 10^(-14), which is
    comparable to the relative precision sqrt implementation of current hardware.
*/
double LPK1(double x);

/** Besselfunctions K0 and K1 at x.
    The implementation has an absolute precision of around 10^(-14), which is
    comparable to the relative precision sqrt implementation of current hardware.
*/
void LPK01(double x, double *K0, double *K1);
#endif
