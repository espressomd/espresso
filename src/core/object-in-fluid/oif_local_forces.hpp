/*
  Copyright (C) 2012-2018 The ESPResSo project

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
#ifndef _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H
#define _OBJECT_IN_FLUID_OIF_LOCAL_FORCES_H

/** \file oif_local_forces.hpp
 *  Routines to calculate the OIF_LOCAL_FORCES
 *  for a particle quadruple (two neighboring triangles with common edge).
 * (Dupin2007) \ref forces.cpp
 */

#include "config.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

// set parameters for local forces
int oif_local_forces_set_params(int bond_type, double r0, double ks,
                                double kslin, double phi0, double kb,
                                double A01, double A02, double kal,
                                double kvisc);

inline double KS(double lambda) { // Defined by (19) from Dupin2007
  double res;
  res = (pow(lambda, 0.5) + pow(lambda, -2.5)) / (lambda + pow(lambda, -3.));
  return res;
}

/** Computes the local forces (Dupin2007) and adds them
    to the particle forces.
    @param p1,p2,p3     Pointers to particles of triangle 1.
    @param p2,p3,p4     Pointers to particles of triangle 2.
    (triangles have particles p2 and p3 in common)
    @return 0
*/
inline int calc_oif_local(Particle *p2, Particle *p1, Particle *p3,
                          Particle *p4, Bonded_ia_parameters *iaparams,
                          double force[3], double force2[3], double force3[3],
                          double force4[3]) // first-fold-then-the-same approach
{
  int i, img[3];
  Vector3d fp1, fp2, fp3, fp4;
  double AA[3], BB[3], CC[3];
  double n1[3], n2[3], dn1, dn2, phi, aa;
  double dx[3], fac, dr, len2, len, lambda;
  double A, h[3], rh[3], hn;
  double m1[3], m2[3], m3[3];
  double v[3], def_vel;
  double m1_length, m2_length, m3_length, t;

  // first find out which particle out of p1, p2 (possibly p3, p4) is not a
  // ghost particle. In almost all cases it is p2, however, it might be other
  // one. we call this particle reference particle.
  if (p2->l.ghost != 1) {
    // unfold non-ghost particle using image, because for physical particles,
    // the structure p->l.i is correctly set
    fp2 = unfolded_position(*p2);
    // other coordinates are obtained from its relative positions to the
    // reference particle
    get_mi_vector(AA, p1->r.p, fp2);
    get_mi_vector(BB, p3->r.p, fp2);
    get_mi_vector(CC, p4->r.p, fp2);
    for (i = 0; i < 3; i++) {
      fp1[i] = fp2[i] + AA[i];
      fp3[i] = fp2[i] + BB[i];
      fp4[i] = fp2[i] + CC[i];
    }
  } else {
    // in case  particle p2 is a ghost particle
    if (p1->l.ghost != 1) {
      fp1 = unfolded_position(*p1);
      get_mi_vector(AA, p2->r.p, fp1);
      get_mi_vector(BB, p3->r.p, fp1);
      get_mi_vector(CC, p4->r.p, fp1);
      for (i = 0; i < 3; i++) {
        fp2[i] = fp1[i] + AA[i];
        fp3[i] = fp1[i] + BB[i];
        fp4[i] = fp1[i] + CC[i];
      }
    } else {
      // in case the first and the second particle are ghost particles
      if (p3->l.ghost != 1) {
        fp3 = unfolded_position(p3);
        get_mi_vector(AA, p1->r.p, fp3);
        get_mi_vector(BB, p2->r.p, fp3);
        get_mi_vector(CC, p4->r.p, fp3);
        for (i = 0; i < 3; i++) {
          fp1[i] = fp3[i] + AA[i];
          fp2[i] = fp3[i] + BB[i];
          fp4[i] = fp3[i] + CC[i];
        }
      } else {
        // in case the first and the second particle are ghost particles
        if (p4->l.ghost != 1) {
          fp4 = unfolded_position(p4);
          get_mi_vector(AA, p1->r.p, fp4);
          get_mi_vector(BB, p2->r.p, fp4);
          get_mi_vector(CC, p3->r.p, fp4);
          for (i = 0; i < 3; i++) {
            fp1[i] = fp4[i] + AA[i];
            fp2[i] = fp4[i] + BB[i];
            fp3[i] = fp4[i] + CC[i];
          }
        } else {
          throw std::runtime_error(
              "Something wrong in oif_local_forces.hpp: All particles in a "
              "bond are ghost "
              "particles, impossible to unfold the positions...\n");
          return 0;
        }
      }
    }
  }

  for (i = 0; i < 3; i++) {
    force[i] = 0;
    force2[i] = 0;
    force3[i] = 0;
    force4[i] = 0;
  }

  // non-linear stretching
  if (iaparams->p.oif_local_forces.ks > TINY_OIF_ELASTICITY_COEFFICIENT) {
    vecsub(fp2, fp3, dx);
    len2 = sqrlen(dx);
    len = sqrt(len2);
    dr = len - iaparams->p.oif_local_forces.r0;
    lambda = 1.0 * len / iaparams->p.oif_local_forces.r0;
    fac =
        -iaparams->p.oif_local_forces.ks * KS(lambda) * dr; // no normalization
    for (i = 0; i < 3; i++) {
      force2[i] += fac * dx[i] / len;
      force3[i] += -fac * dx[i] / len;
    }
  }

  // linear stretching
  if (iaparams->p.oif_local_forces.kslin > TINY_OIF_ELASTICITY_COEFFICIENT) {
    vecsub(fp2, fp3, dx);
    len2 = sqrlen(dx);
    len = sqrt(len2);
    dr = len - iaparams->p.oif_local_forces.r0;
    fac = -iaparams->p.oif_local_forces.kslin * dr; // no normalization
    for (i = 0; i < 3; i++) {
      force2[i] += fac * dx[i] / len;
      force3[i] += -fac * dx[i] / len;
    }
  }

  // viscous force
  if (iaparams->p.oif_local_forces.kvisc >
      TINY_OIF_ELASTICITY_COEFFICIENT) { // to be implemented....

    vecsub(fp2, fp3, dx);
    len2 = sqrlen(dx);
    len = sqrt(len2);

    v[0] = (p3->m.v[0] - p2->m.v[0]);
    v[1] = (p3->m.v[1] - p2->m.v[1]);
    v[2] = (p3->m.v[2] - p2->m.v[2]);

    // Variant A
    // Here the force is in the direction of relative velocity btw points

    // Code:
    // for(i=0;i<3;i++) {
    // force2[i] += iaparams->p.oif_local_forces.kvisc*v[i];
    // force3[i] -= iaparams->p.oif_local_forces.kvisc*v[i];
    //}

    // Variant B
    // Here the force is the projection of relative velocity btw points onto
    // line btw the points

    // denote p vector between p2 and p3
    // denote v the velocity difference between the points p2 and p3
    // denote alpha the angle between p and v
    // denote x the projected v onto p
    // cos alpha = |x|/|v|
    // cos alpha = (v,p)/(|v||p|)
    // together we get |x|=(v,p)/|p|
    // also, x is along p, so x = |x|.p/|p|
    // so x = p/|p| . (v,p)/|p|
    // altogether x = p . (v,p)/(|p|)^2
    // |p|^2 is stored in len2

    // Code:
    def_vel = dx[0] * v[0] + dx[1] * v[1] + dx[2] * v[2];
    fac = iaparams->p.oif_local_forces.kvisc * def_vel / len2;
    for (i = 0; i < 3; i++) {
      force2[i] += fac * dx[i];
      force3[i] -= fac * dx[i];
    }
  }

  /* bending
     forceT1 is restoring force for triangle p1,p2,p3 and force2T restoring
     force for triangle p2,p3,p4 p1 += forceT1; p2 -= 0.5*forceT1+0.5*forceT2;
     p3 -= 0.5*forceT1+0.5*forceT2; p4 += forceT2; */
  if (iaparams->p.oif_local_forces.kb > TINY_OIF_ELASTICITY_COEFFICIENT) {
    get_n_triangle(fp2, fp1, fp3, n1);
    dn1 = normr(n1);
    get_n_triangle(fp2, fp3, fp4, n2);
    dn2 = normr(n2);
    phi = angle_btw_triangles(fp1, fp2, fp3, fp4);

    aa = (phi - iaparams->p.oif_local_forces
                    .phi0); // no renormalization by phi0, to be consistent with
                            // Krueger and Fedosov
    fac = iaparams->p.oif_local_forces.kb * aa;
    for (i = 0; i < 3; i++) {
      force[i] += fac * n1[i] / dn1;
      force2[i] -= (0.5 * fac * n1[i] / dn1 + 0.5 * fac * n2[i] / dn2);
      force3[i] -= (0.5 * fac * n1[i] / dn1 + 0.5 * fac * n2[i] / dn2);
      force4[i] += fac * n2[i] / dn2;
    }
  }

  /* local area
     for both triangles
     only 1/3 of calculated forces are added, because each triangle will enter
     this calculation 3 times (one time per edge)

              Proportional distribution of forces, implemented according to the
     article I.Jancigova, I.Cimrak, Non-uniform force allocation for area
     preservation in spring network models, International Journal for Numerical
     Methods in Biomedical Engineering, DOI: 10.1002/cnm.2757

  */
  if (iaparams->p.oif_local_forces.kal > TINY_OIF_ELASTICITY_COEFFICIENT) {

    // triangle p1,p2,p3
    for (i = 0; i < 3; i++) { // centroid of triangle p1,p2,p3
      h[i] = 1.0 / 3.0 * (fp1[i] + fp2[i] + fp3[i]);
    }
    A = area_triangle(fp1, fp2, fp3);
    t = sqrt(A / iaparams->p.oif_local_forces.A01) - 1.0;
    vecsub(h, fp1, m1);
    vecsub(h, fp2, m2);
    vecsub(h, fp3, m3);

    m1_length = normr(m1);
    m2_length = normr(m2);
    m3_length = normr(m3);

    fac =
        iaparams->p.oif_local_forces.kal * iaparams->p.oif_local_forces.A01 *
        (2 * t + t * t) /
        (m1_length * m1_length + m2_length * m2_length + m3_length * m3_length);

    for (i = 0; i < 3; i++) { // local area force for p1
      force[i] += fac * m1[i] / 3.0;
    }
    for (i = 0; i < 3; i++) { // local area force for p2
      force2[i] += fac * m2[i] / 3.0;
    }
    for (i = 0; i < 3; i++) { // local area force for p3
      force3[i] += fac * m3[i] / 3.0;
    }

    // triangle p2,p3,p4
    for (i = 0; i < 3; i++) { // centroid of triangle p2,p3,p4
      h[i] = 1.0 / 3.0 * (fp2[i] + fp3[i] + fp4[i]);
    }
    A = area_triangle(fp2, fp3, fp4);
    t = sqrt(A / iaparams->p.oif_local_forces.A02) - 1.0; ////
    vecsub(h, fp2, m1);
    vecsub(h, fp3, m2);
    vecsub(h, fp4, m3);

    m1_length = normr(m1);
    m2_length = normr(m2);
    m3_length = normr(m3);

    fac =
        iaparams->p.oif_local_forces.kal * iaparams->p.oif_local_forces.A02 *
        (2 * t + t * t) /
        (m1_length * m1_length + m2_length * m2_length + m3_length * m3_length);

    for (i = 0; i < 3; i++) { // local area force for p2
      force2[i] += fac * m1[i] / 3.0;
    }
    for (i = 0; i < 3; i++) { // local area force for p3
      force3[i] += fac * m2[i] / 3.0;
    }
    for (i = 0; i < 3; i++) { // local area force for p4
      force4[i] += fac * m3[i] / 3.0;
    }
  }
  return 0;
}

#endif
