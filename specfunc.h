// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef SPECFUNC_H
#define SPECFUNC_H

/** Hurwitz zeta function. This function was taken from the GSL code. */
double hzeta(double order, double x);

/** Modified Bessel function of first kind, order 0. This function was taken from
    the GSL code. */
double I0(double x);
/** Modified Bessel function of first kind, order 1. This function was taken from
    the GSL code. */
double I1(double x);
/** Modified Bessel function of second kind, order 0. This function was taken from
    the GSL code. */
double K0(double x);
/** Modified Bessel function of second kind, order 1. This function was taken from
    the GSL code. */
double K1(double x);

#endif
