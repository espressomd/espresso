/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "Ellipsoid.hpp"
#include "errorhandling.hpp"
#include "utils.hpp"
#include <cmath>

namespace Shapes {
int Ellipsoid::calculate_dist(const double *ppos, double *dist,
                              double *vec) const {

  /* get particle position in reference frame of ellipsoid */
  Vector3d const ppos_e = Vector3d(ppos, ppos + 3) - m_center;

  /* set appropriate initial point for Newton's method */
  double l0, l = 0.;
  int distance_prefactor = -1;
  if (not inside_ellipsoid(ppos_e)) {
    l = *std::max_element(m_semiaxes.begin(), m_semiaxes.end()) * ppos_e.norm();
    distance_prefactor = 1;
  }

  /* find root via Newton's method */
  double eps = 10.;
  int step = 0;
  while ((eps >= 1e-12) and (step < 100)) {
    l0 = l;
    l -= newton_term(ppos_e, l0);
    eps = std::abs(l - l0);
    step++;
  }

  if (step == 100) {
    runtimeWarning("Maximum number of Newton steps exceeded while calculating "
                   "distance to Ellipsoid.");
  }

  /* calculate dist and vec */
  double distance = 0.;
  for (int i = 0; i < 3; i++) {
    vec[i] = (ppos_e[i] -
              Utils::sqr(m_semiaxes[i]) * ppos_e[i] /
                  (l + Utils::sqr(m_semiaxes[i])));
    distance += Utils::sqr(vec[i]);
  }

  *dist = distance_prefactor * m_direction * std::sqrt(distance);

  return 0;
}

bool Ellipsoid::inside_ellipsoid(const Vector3d &ppos) const {
  bool is_inside = false;
  if (Utils::sqr(ppos[0] / m_semiaxes[0]) +
          Utils::sqr(ppos[1] / m_semiaxes[1]) +
          Utils::sqr(ppos[2] / m_semiaxes[2]) <=
      1)
    is_inside = true;

  return is_inside;
}

double Ellipsoid::newton_term(const Vector3d &ppos, const double &l) const {
  Vector3d axpos, lax, lax2;
  for (int i = 0; i < 3; i++) {
    axpos[i] = Utils::sqr(m_semiaxes[i]) * Utils::sqr(ppos[i]);
    lax[i] = l + Utils::sqr(m_semiaxes[i]);
    lax2[i] = Utils::sqr(lax[i]);
  }

  return (axpos[0] * lax2[1] * lax2[2] + axpos[1] * lax2[2] * lax2[0] +
          axpos[2] * lax2[0] * lax2[1] - lax2[0] * lax2[1] * lax2[2]) /
         (2 * (axpos[0] * lax[1] * lax2[2] + axpos[1] * lax[2] * lax2[0] +
               axpos[2] * lax[0] * lax2[1] + axpos[0] * lax2[1] * lax[2] +
               axpos[1] * lax2[2] * lax[0] + axpos[2] * lax2[0] * lax[1] -
               lax[0] * lax2[1] * lax2[2] - lax2[0] * lax[1] * lax2[2] -
               lax2[0] * lax2[1] * lax[2]));
}

} // namespace Shapes
