#ifndef UTILS_MATH_TRIANGEL_FUNCTIONS_HPP
#define UTILS_MATH_TRIANGEL_FUNCTIONS_HPP

#include "utils/Vector.hpp"

#include <cmath>

namespace Utils {
/** Computes the normal vector to the plane given by points P1P2P3 */
inline Vector3d get_n_triangle(const Vector3d &P1, const Vector3d &P2,
                               const Vector3d &P3) {
  auto const u = P2 - P1;
  auto const v = P3 - P1;

  return vector_product(u, v);
}

/** Computes the area of triangle between vectors P1,P2,P3,
 *  by computing the crossproduct P1P2 x P1P3 and taking the half of its norm */
inline double area_triangle(const Vector3d &P1, const Vector3d &P2,
                            const Vector3d &P3) {
  return 0.5 * get_n_triangle(P1, P2, P3).norm();
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
double angle_btw_triangles(const T1 &P1, const T2 &P2, const T3 &P3,
                           const T4 &P4) {
  auto const normal1 = get_n_triangle(P2, P1, P3);
  auto const normal2 = get_n_triangle(P2, P3, P4);

  double tmp11;
  // Now we compute the scalar product of n1 and n2 divided by the norms of n1
  // and n2
  // tmp11 = dot(3,normal1,normal2);         // tmp11 = n1.n2
  tmp11 = normal1 * normal2; // tmp11 = n1.n2

  tmp11 *= fabs(tmp11);                         // tmp11 = (n1.n2)^2
  tmp11 /= (normal1.norm2() * normal2.norm2()); // tmp11 = (n1.n2/(|n1||n2|))^2
  if (tmp11 > 0) {
    tmp11 = std::sqrt(tmp11);
  } else {
    tmp11 = -std::sqrt(-tmp11);
  }

  if (tmp11 >= 1.) {
    tmp11 = 0.0;
  } else if (tmp11 <= -1.) {
    tmp11 = M_PI;
  }
  auto const phi =
      M_PI - std::acos(tmp11); // The angle between the faces (not considering
                               // the orientation, always less or equal to Pi)
                               // is equal to Pi minus angle between the normals

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
    return 2 * M_PI - phi;

  return phi;
}
} // namespace Utils

#endif
