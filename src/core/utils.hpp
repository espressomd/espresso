/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef UTILS_HPP
#define UTILS_HPP
/** \file utils.hpp
 *    Small functions that are useful not only for one modul.

 *  just some nice utilities...
 *  and some constants...
 *
*/

#include "Vector.hpp"
#include "utils/constants.hpp"
#include "utils/math/sqr.hpp"
#include "utils/memory.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

namespace Utils {
/**
 * \brief Calculate integer powers.
 * This functions calculates x^n, where
 * n is a positive integer that is known
 * at compile time. It uses exponentiation by
 * squaring to construct a efficient function.
 */
template <unsigned n, typename T> inline T int_pow(T x) {
  switch (n) {
  case 0:
    return T(1);
  case 1:
    return x;
  default:
    /** Even branch */
    if (n % 2 == 0) {
      return int_pow<n / 2, T>(x * x);
    } else {
      return x * int_pow<(n - 1) / 2, T>(x * x);
    }
  }
}

/** Calculate signum of val, if supported by T */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

} // Namespace Utils

/*************************************************************/
/** \name Mathematical functions.                            */
/*************************************************************/
/*@{*/

/** permute an integer array field of size size about permute positions. */
inline void permute_ifield(int *field, int size, int permute) {
   if (permute == 0)
    return;
  if (permute < 0)
    permute = (size + permute);
  while (permute > 0) {
    int tmp = field[0];
    for (int i = 1; i < size; i++)
      field[i - 1] = field[i];
    field[size - 1] = tmp;
    permute--;
  }
}

/** Mathematically rounds 'double'-typed x, returning 'double'. */
inline double dround(double x) { return floor(x + 0.5); }

/** approximates \f$ \exp(d^2) \mathrm{erfc}(d)\f$ by applying a formula from:
    Abramowitz/Stegun: Handbook of Mathematical Functions, Dover
    (9. ed.), chapter 7 */
inline double AS_erfc_part(double d) {
#define AS_a1 0.254829592
#define AS_a2 -0.284496736
#define AS_a3 1.421413741
#define AS_a4 -1.453152027
#define AS_a5 1.061405429
#define AS_p 0.3275911
  double t;

  t = 1.0 / (1.0 + AS_p * d);

  return t * (AS_a1 + t * (AS_a2 + t * (AS_a3 + t * (AS_a4 + t * AS_a5))));
}

/** Calculates the sinc-function as sin(PI*x)/(PI*x).
 *
 * (same convention as in Hockney/Eastwood). In order to avoid
 * divisions by 0, arguments, whose modulus is smaller than epsi, will
 * be evaluated by an 8th order Taylor expansion of the sinc
 * function. Note that the difference between sinc(x) and this
 * expansion is smaller than 0.235e-12, if x is smaller than 0.1. (The
 * next term in the expansion is the 10th order contribution
 * PI^10/39916800 * x^10 = 0.2346...*x^12).  This expansion should
 * also save time, since it reduces the number of function calls to
 * sin().
 */
inline double sinc(double d) {
  constexpr double epsi = 0.1;

  const double PId = PI * d;

  if (fabs(d) > epsi)
    return sin(PId) / PId;
  else {
    /** Coefficients of the Taylor expansion of sinc */
    constexpr double c2 = -0.1666666666667e-0;
    constexpr double c4 = 0.8333333333333e-2;
    constexpr double c6 = -0.1984126984127e-3;
    constexpr double c8 = 0.2755731922399e-5;

    const double PId2 = PId * PId;
    return 1.0 + PId2 * (c2 + PId2 * (c4 + PId2 * (c6 + PId2 * c8)));
  }
}
/*@}*/

/*************************************************************/
/** \name Vector and matrix operations for three dimensons.  */
/*************************************************************/
/*@{*/

/** Subtracts vector v2 from vector v1 and stores result in vector dv */
template <typename T, typename U, typename V>
inline void vecsub(T const &v1, U const &v2, V &dv) {
  dv[0] = v1[0] - v2[0];
  dv[1] = v1[1] - v2[1];
  dv[2] = v1[2] - v2[2];
}

/** calculates the length of a vector */
inline double normr(double v[3]) {
  double d2 = 0.0;
  int i;
  for (i = 0; i < 3; i++)
    d2 += Utils::sqr(v[i]);
  d2 = sqrtf(d2);
  return d2;
}

/** calculates the squared length of a vector */
template <typename T> double sqrlen(T const &v) {
  double d2 = 0.0;
  int i;
  for (i = 0; i < 3; i++)
    d2 += Utils::sqr(v[i]);
  return d2;
}

/** calculates unit vector */
inline void unit_vector(double v[3], double y[3]) {
  double d = 0.0;
  int i;
  d = sqrt(sqrlen(v));

  for (i = 0; i < 3; i++)
    y[i] = v[i] / d;

  return;
}

/** calculates the scalar product of two vectors a nd b */
template <typename T1, typename T2>
double scalar(const T1& a, const T2& b) {
  double d2 = 0.0;
  int i;
  for (i = 0; i < 3; i++)
    d2 += a[i] * b[i];
  return d2;
}

/** calculates the vector product c of two vectors a and b */
template <typename T, typename U, typename V>
inline void vector_product(T const &a, U const &b, V &c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return;
}

/** rotates vector around axis by alpha */
template <typename T1, typename T2, typename T3>
void vec_rotate(const T1& axis, double alpha, const T2& vector,
                       T3& result) {
  double sina, cosa, absa, a[3];
  sina = sin(alpha);
  cosa = cos(alpha);
  absa = sqrt(scalar(axis, axis));

  a[0] = axis[0] / absa;
  a[1] = axis[1] / absa;
  a[2] = axis[2] / absa;

  result[0] = (cosa + Utils::sqr(a[0]) * (1 - cosa)) * vector[0] +
              (a[0] * a[1] * (1 - cosa) - a[2] * sina) * vector[1] +
              (a[0] * a[2] * (1 - cosa) + a[1] * sina) * vector[2];
  result[1] = (a[0] * a[1] * (1 - cosa) + a[2] * sina) * vector[0] +
              (cosa + Utils::sqr(a[1]) * (1 - cosa)) * vector[1] +
              (a[1] * a[2] * (1 - cosa) - a[0] * sina) * vector[2];
  result[2] = (a[0] * a[2] * (1 - cosa) - a[1] * sina) * vector[0] +
              (a[1] * a[2] * (1 - cosa) + a[0] * sina) * vector[1] +
              (cosa + Utils::sqr(a[2]) * (1 - cosa)) * vector[2];

  return;
}

/** rotates vector around axis by alpha */
inline ::Vector<3, double> vec_rotate(::Vector<3, double> axis, double alpha,
                                      ::Vector<3, double> vector) {
  ::Vector<3, double> result;
  vec_rotate(axis, alpha, vector, result);
  return result;
}

/*@}*/

/*************************************************************/
/** \name Three dimensional grid operations                  */
/*************************************************************/
/*@{*/

/** get the linear index from the position (a,b,c) in a 3D grid
 *  of dimensions adim[]. returns linear index.
 *
 * @return        the linear index
 * @param a       x position
 * @param b       y position
 * @param c       z position
 * @param adim    dimensions of the underlying grid
 */
inline int get_linear_index(int a, int b, int c, int adim[3]) {
  return (a + adim[0] * (b + adim[1] * c));
}

/** get the position a[] from the linear index in a 3D grid
 *  of dimensions adim[].
 *
 * @param i       linear index
 * @param a       x position (return value)
 * @param b       y position (return value)
 * @param c       z position (return value)
 * @param adim    dimensions of the underlying grid
 */
inline void get_grid_pos(int i, int *a, int *b, int *c, int adim[3]) {
  *a = i % adim[0];
  i /= adim[0];
  *b = i % adim[1];
  i /= adim[1];
  *c = i;
}

/*@}*/

/*************************************************************/
/** \name Distance calculations.  */
/*************************************************************/
/*@{*/

/** returns the distance between two position.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 */
inline double distance(double pos1[3], double pos2[3]) {
  return sqrt(Utils::sqr(pos1[0] - pos2[0]) + Utils::sqr(pos1[1] - pos2[1]) +
              Utils::sqr(pos1[2] - pos2[2]));
}

/** returns the distance between two positions squared.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 */
template <typename T1, typename T2>
inline double distance2(const T1& pos1, const T2& pos2) {
  return Utils::sqr(pos1[0] - pos2[0]) + Utils::sqr(pos1[1] - pos2[1]) +
         Utils::sqr(pos1[2] - pos2[2]);
}

/** Returns the distance between two positions squared and stores the
    distance vector pos1-pos2 in vec.
 *  \param pos1 Position one.
 *  \param pos2 Position two.
 *  \param vec  vecotr pos1-pos2.
 *  \return distance squared
*/
template <typename T1, typename T2, typename T3>
double distance2vec(T1 const pos1, T2 const pos2,
                           T3& vec) {
  vec[0] = pos1[0] - pos2[0];
  vec[1] = pos1[1] - pos2[1];
  vec[2] = pos1[2] - pos2[2];
  return Utils::sqr(vec[0]) + Utils::sqr(vec[1]) + Utils::sqr(vec[2]);
}

/*@}*/

/*************************************************************/
/** \name String helper functions                            */
/*************************************************************/
/*@{*/

/** extend a string with another one. Like strcat, just automatically
    increases the string space */
char *strcat_alloc(char *left, const char *right);

/*@}*/

/*************************************************************/
/** \name Object-in-fluid functions                          */
/*************************************************************/
/*@{*/

/** Computes the area of triangle between vectors P1,P2,P3,
 *  by computing the crossproduct P1P2 x P1P3 and taking the half of its norm */
template <typename T1, typename T2, typename T3>
inline double area_triangle(const T1& P1, const T2& P2, const T3& P3) {
  double area;
  double u[3], v[3], normal[3], n; // auxiliary variables
  u[0] = P2[0] - P1[0];            // u = P1P2
  u[1] = P2[1] - P1[1];
  u[2] = P2[2] - P1[2];
  v[0] = P3[0] - P1[0]; // v = P1P3
  v[1] = P3[1] - P1[1];
  v[2] = P3[2] - P1[2];
  vector_product(u, v, normal);
  n = normr(normal);
  area = 0.5 * n;
  return area;
}

/** Computes the normal vector to the plane given by points P1P2P3 */
template <typename T1, typename T2, typename T3>
void get_n_triangle(const T1& p1, const T2& p2, const T3& p3, double *n) {
  n[0] = (p2[1] - p1[1]) * (p3[2] - p1[2]) - (p2[2] - p1[2]) * (p3[1] - p1[1]);
  n[1] = (p2[2] - p1[2]) * (p3[0] - p1[0]) - (p2[0] - p1[0]) * (p3[2] - p1[2]);
  n[2] = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
}

/** This function returns the angle btw the triangle p1,p2,p3 and p2,p3,p4.  Be
 * careful, the angle depends on the orientation of the trianlges!  You need to
 * be sure that the orientation (direction of normal vector) of p1p2p3 is given
 * by the cross product p2p1 x p2p3.  The orientation of p2p3p4 must be given
 * by p2p3 x p2p4.
 *
 *  Example: p1 = (0,0,1), p2 = (0,0,0), p3=(1,0,0), p4=(0,1,0).  The
 *  orientation of p1p2p3 should be in the direction (0,1,0) and indeed: p2p1 x
 *  p2p3 = (0,0,1)x(1,0,0) = (0,1,0) This function is called in the beginning
 *  of the simulation when creating bonds depending on the angle btw the
 *  triangles, the bending_force.  Here, we determine the orientations by
 *  looping over the triangles and checking the correct orientation.  So if you
 *  have the access to the order of particles, you are safe to call this
 *  function with exactly this order. Otherwise you need to check the
 *  orientations. */
template <typename T1, typename T2, typename T3, typename T4>
double angle_btw_triangles(const T1& P1, const T2& P2, const T3& P3, const T4& P4) {
  double phi;
  double u[3], v[3], normal1[3], normal2[3]; // auxiliary variables
  u[0] = P1[0] - P2[0];                      // u = P2P1
  u[1] = P1[1] - P2[1];
  u[2] = P1[2] - P2[2];
  v[0] = P3[0] - P2[0]; // v = P2P3
  v[1] = P3[1] - P2[1];
  v[2] = P3[2] - P2[2];
  vector_product(u, v, normal1);
  u[0] = P3[0] - P2[0]; // u = P2P3
  u[1] = P3[1] - P2[1];
  u[2] = P3[2] - P2[2];
  v[0] = P4[0] - P2[0]; // v = P2P4
  v[1] = P4[1] - P2[1];
  v[2] = P4[2] - P2[2];
  vector_product(u, v, normal2);

  double tmp11, tmp22, tmp33;
  // Now we compute the scalar product of n1 and n2 divided by the norms of n1
  // and n2
  // tmp11 = dot(3,normal1,normal2);         // tmp11 = n1.n2
  tmp11 = scalar(normal1, normal2); // tmp11 = n1.n2

  /*
  tmp22 = normr(normal1);
  tmp33 = normr(normal2);
  tmp11 /= (tmp22*tmp33);  // tmp11 = n1.n2/(|n1||n2|)
*/
  tmp11 *= fabs(tmp11);     // tmp11 = (n1.n2)^2
  tmp22 = sqrlen(normal1);  // tmp22 = |n1|^2
  tmp33 = sqrlen(normal2);  // tmp33 = |n1|^2
  tmp11 /= (tmp22 * tmp33); // tmp11 = (n1.n2/(|n1||n2|))^2
  if (tmp11 > 0) {
    tmp11 = sqrt(tmp11);
  } else {
    tmp11 = -sqrt(-tmp11);
  }

  if (tmp11 >= 1.) {
    tmp11 = 0.0;
  } else if (tmp11 <= -1.) {
    tmp11 = M_PI;
  }
  phi = M_PI - acos(tmp11); // The angle between the faces (not considering the
                            // orientation, always less or equal to Pi) is
                            // equal to Pi minus angle between the normals

  // Now we need to determine, if the angle btw two triangles is less than Pi or
  // more than Pi. To do this we check,
  // if the point P4 lies in the halfspace given by trianlge P1P2P3 and the
  // normal to this triangle. If yes, we have
  // angle less than Pi, if not, we have angle more than Pi.
  // General equation of the plane is n_x*x + n_y*y + n_z*z + d = 0 where
  // (n_x,n_y,n_z) is the normal to the plane.
  // Point P1 lies in the plane, therefore d = -(n_x*P1_x + n_y*P1_y + n_z*P1_z)
  // Point P4 lies in the halfspace given by normal iff n_x*P4_x + n_y*P4_y +
  // n_z*P4_z + d >= 0
  tmp11 = -(normal1[0] * P1[0] + normal1[1] * P1[1] + normal1[2] * P1[2]);
  if (normal1[0] * P4[0] + normal1[1] * P4[1] + normal1[2] * P4[2] + tmp11 < 0)
    phi = 2 * M_PI - phi;
  return phi;
}

namespace Utils {

// Below you will find some routines for handling vectors.  Note that
// all the pointer-type overloads assume a pointer of length 3!  Due
// to restrictions on the C-level there is no error checking available
// for the pointer-type overloads, i.e. you won't trigger an exception
// in out-of-memory situations but will receive a segfault right away.

//
// cross_product: Calculate the cross product of two vectors
//
template <typename T1, typename T2, typename T3>
void cross_product(const T1& a, const T2& b, T3& c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

template <typename T> double dot_product(T const *const a, T const *const b) {
  double c = 0;
  for (unsigned int i = 0; i < 3; i++)
    c += a[i] * b[i];
  return c;
}

//
// veclen and sqrlen: Calculate the length and length squared of a vector
//
template <typename T> double sqrlen(T const *const a) {
  double c = 0;
  for (int i = 0; i < 3; ++i)
    c += a[i] * a[i];
  return c;
}

template <typename T> double veclen(T const *const a) {
  return sqrt(sqrlen(a));
}

//
// vecsub: Subtract two vectors
//

template <typename T>
void vecsub(T const *const a, T const *const b, T *const c) {
  // Note the different signature for pointers here!
  for (unsigned int i = 0; i < 3; i++)
    c[i] = a[i] - b[i];
}

} // namespace Utils

/*@}*/

#endif
