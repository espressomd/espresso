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
