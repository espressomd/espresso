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
/** \file angledist.cpp
 *
 *  Implementation of \ref angledist.hpp
 */
#include "angledist.hpp"

#ifdef BOND_ANGLEDIST

#include <algorithm>

#include "communication.hpp"
#include "constraints.hpp"
#include "constraints/ShapeBasedConstraint.hpp"
#include "grid.hpp"
#include <memory>
#include "shapes/Wall.hpp"

int angledist_set_params(int bond_type, double bend, double phimin,
                         double distmin, double phimax, double distmax) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angledist.bend = bend;
  bonded_ia_params[bond_type].p.angledist.phimin = phimin;
  bonded_ia_params[bond_type].p.angledist.distmin = distmin;
  bonded_ia_params[bond_type].p.angledist.phimax = phimax;
  bonded_ia_params[bond_type].p.angledist.distmax = distmax;
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLEDIST;
  bonded_ia_params[bond_type].num = 2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

/** Function to calculate wall distance and phi0(dist).
    Called by \ref calc_angledist_force */
static double calc_angledist_param(Particle *p_mid, Particle *p_left,
                                   Particle *p_right,
                                   Bonded_ia_parameters *iaparams) {
  double vec1[3], vec2[3], d2i = 0.0, dist2 = 0.0;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  const double dist1 = sqrlen(vec1);
  const double d1i = 1.0 / sqrt(dist1);
  for (int j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (int j = 0; j < 3; j++)
    vec2[j] *= d2i;

  const double phimn = iaparams->p.angledist.phimin;
  const double distmn = iaparams->p.angledist.distmin;
  const double phimx = iaparams->p.angledist.phimax;
  const double distmx = iaparams->p.angledist.distmax;

  /* folds coordinates of p_mid into original box */
  Vector3d folded_pos =folded_position(*p_mid);

  double pwdistmin = std::numeric_limits<double>::infinity();
  double pwdistmin_d = 0.0;

  for (auto const &c : Constraints::constraints) {
    auto cs = std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
    if (cs) {
      if (dynamic_cast<const Shapes::Wall*>(&(cs->shape()))) {
        const Shapes::Wall* wall = dynamic_cast<const Shapes::Wall*>(&(cs->shape()));
        double dist = -wall->d();
        for (int j = 0; j < 3; j++) {
          dist += folded_pos[j] * (wall->n())[j];
        }

        if (dist < pwdistmin) {
          pwdistmin = dist;
          pwdistmin_d = wall->d();
        }
      }
    }
  }

  /*get phi0(z)*/
  if (pwdistmin <= distmn) {
    return phimn;
  } else if (pwdistmin >= distmx &&
             pwdistmin <= box_l[2] - pwdistmin_d - distmx) {
    return phimx;
  } else {
    const double drange = (pwdistmin - distmn) * PI / (distmx - distmn);
    return ((cos(drange - PI) + 1.0) * (phimx - phimn)) * 0.5 + phimn;
  }
}

int calc_angledist_force(Particle *p_mid, Particle *p_left, Particle *p_right,
                         Bonded_ia_parameters *iaparams, double force1[3],
                         double force2[3]) {
  double cosine = 0.0, vec1[3], vec2[3], fac = 0.0;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  const double dist12 = sqrlen(vec1);
  const double d1i = 1.0 / sqrt(dist12);
  for (int j = 0; j < 3; j++)
    vec1[j] *= d1i;

  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  const double dist22 = sqrlen(vec2);
  const double d2i = 1.0 / sqrt(dist22);
  for (int j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar product of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  fac = iaparams->p.angledist.bend;

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    const double phi0 = calc_angledist_param(p_mid, p_left, p_right, iaparams);

    if (cosine > TINY_COS_VALUE)
      cosine = TINY_COS_VALUE;
    if (cosine < -TINY_COS_VALUE)
      cosine = -TINY_COS_VALUE;
    const double phi = acos(-cosine);

    double sinphi = sin(phi);
    if (sinphi < TINY_SIN_VALUE)
      sinphi = TINY_SIN_VALUE;
    fac *= (phi - phi0) / sinphi;
  }
#endif

  for (int j = 0; j < 3; j++) {
    double f1 = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    double f2 = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force1[j] = (f1 - f2);
    force2[j] = -f1;
  }
  return 0;
}

#ifdef BOND_ANGLEDIST_HARMONIC
int angledist_energy(Particle *p_mid, Particle *p_left, Particle *p_right,
                     Bonded_ia_parameters *iaparams, double *_energy) {
  int j;
  double cosine = 0.0, d1i, d2i, dist1, dist2;
  double vec1[3], vec2[3];

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for (j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;

  {
    double phi;
    double phi0 = calc_angledist_param(p_mid, p_left, p_right, iaparams);
    phi = acos(-cosine);
    *_energy = 0.5 * iaparams->p.angledist.bend * Utils::sqr(phi - phi0);
  }

  return 0;
}
#endif

#endif
