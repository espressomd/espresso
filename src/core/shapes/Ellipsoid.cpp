/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
#include <cmath>
#include <iostream>

using namespace std;

#define SQR(A) ((A) * (A))

namespace Shapes {
int Ellipsoid::calculate_dist(const double *ppos, double *dist,
                           double *vec) const {

    Vector3d pposC, pposQ1, pelQ1, pel;

    // Shift ppos into reference frame of ellipsoid
    for (int i=0; i<3; i++) {
        pposC[i] = ppos[i] - m_pos[i];
    }

    // Flip into first Quadrant
    for (int i=0; i<3; i++) {
        pposQ1[i] = abs(pposC[i]);
    }

    // Calc closest point
    pelQ1 = ClosestEllipsoidPoint(pposQ1);

    // Flip back if necessary
    for (int i=0; i<3; i++) {
        pel[i] = std::copysign(pelQ1[i], pposC[i]);
    }

    double distance = 0;

    for (int i=0; i<3; i++) {
        vec[i] = m_direction * (pel[i] - pposC[i]);
        distance += SQR(pel[i] - pposC[i]);
    } 

    *dist = sqrt(distance);

    return 0;
}

// Requires that a, b, c > 0 and ppos[i] > 0 for all i
Vector3d Ellipsoid::ClosestEllipsoidPoint(Vector3d ppos) const {

    Vector3d closest_point;

    if (m_semiaxis_c > 0) {
        if (m_semiaxis_b > 0) {
            if (m_semiaxis_a > 0) {
                const double za = ppos[0]/m_semiaxis_a;
                const double zb = ppos[1]/m_semiaxis_b;
                const double zc = ppos[2]/m_semiaxis_c;
                const double g = SQR(za) + SQR(zb) + SQR(zc) - 1;
                if (g != 0) {
                    const double ra = SQR(m_semiaxis_a/m_semiaxis_c);
                    const double rb = SQR(m_semiaxis_b/m_semiaxis_c);
                    double sbar = GetRoot(ra, rb, za, zb, zc, g);
                    closest_point[0] = ra * ppos[0]/(sbar + ra);
                    closest_point[1] = rb * ppos[1]/(sbar + rb);
                    closest_point[2] = ppos[2]/(sbar + 1);
                }
                else {
                    for (int i=0; i<3; i++) {
                        closest_point[i] = ppos[i];
                    }
                }
            }
            else { // ya == 0
                closest_point[0] = 0;
            }
        }
        else { // yb == 0
            if (ppos[0] > 0) {
                closest_point[1] = 0;
            }
            else {
                closest_point[0] = 0;
                closest_point[1] = 0;
                closest_point[2] = m_semiaxis_c;
            }
        }
    }
    else { // yc == 0
        double denoma = m_semiaxis_a * m_semiaxis_a - m_semiaxis_c * m_semiaxis_c;
        double denomb = m_semiaxis_b * m_semiaxis_b - m_semiaxis_c * m_semiaxis_c;
        double numera = m_semiaxis_a * ppos[0];
        double numerb = m_semiaxis_b * ppos[1];
        bool computed = false;
        if (numera < denoma && numerb < denomb) {
            double xdea = numera/denoma;
            double xdeb = numerb/denomb;
            double xdeasqr = xdea * xdea;
            double xdebsqr = xdeb * xdeb;
            double discr = 1 - xdeasqr - xdebsqr;
            if (discr > 0) {
                closest_point[0] = m_semiaxis_a * xdea;
                closest_point[1] = m_semiaxis_b * xdeb;
                closest_point[2] = m_semiaxis_c * sqrt(discr);
                computed = true;
            }
        }
        if (not computed) {
            closest_point[2] = 0;
        }
    }
    return closest_point;

}

double Ellipsoid::DistancePointEllipse(const double e0, const double e1, const double y0, const double y1, double& x0, double& x1) const {
    double distance;
    if (y1 > 0) {
        if (y0 > 0) {
            double z0 = y0/e0;
            double z1 = y1/e1;
            double g = z0*z0 + z1*z1 - 1;
            if (g != 0) {
                double r0 = SQR(e0/e1);
                double sbar = GetRoot(r0, z0, z1, g);
                x0 = r0*y0/(sbar + r0);
                x1 = y1/(sbar + 1);
                distance = sqrt(SQR(x0 - y0) + SQR(x1 - y1));
            }
            else {
                x0 = y0;
                x1 = y1;
                distance = 0;
            }
        }
        else { // y0 == 0
            x0 = 0;
            x1 = e1;
            distance = abs (y1 - e1);
        }
    }
    else { // y1 == 0
        double numer0 = e0*y0;
        double denom0 = SQR(e0) - SQR(e1);
        if (numer0 < denom0) {
            double xde0 = numer0/denom0;
            x0 = e0*xde0;
            x1 = e1*sqrt(1-xde0*xde0);
            distance = sqrt(SQR(x0-y0) + SQR(x1));
        }
        else {
            x0 = e0;
            x1 = 0;
            distance = abs(y0 - e0);
        }
    }
    return distance;
}

// For Ellipsoid
double Ellipsoid::RobustLength(const double v0, const double v1, const double v2) const {
    double vmax = std::max(v0, std::max(v1, v2));
    return abs(vmax)*sqrt(SQR(v0/vmax) + SQR(v1/vmax) + SQR(v2/vmax));
}

// For Ellipsoid
double Ellipsoid::GetRoot(const double r0, const double r1, const double z0, const double z1, const double z2, double g) const {
    double n0 = r0*z0;
    double n1 = r1*z1;
    double s0 = z1 - 1;
    double s1 = (g < 0 ? 0 : RobustLength(n0, n1, z2) - 1);
    double s = 0;
    int maxIterations = 100;
    for (int i=0; i<maxIterations; ++i) {
        s = (s0+s1)/2;
        if (s == s0 or s == s1) { break; }
        double ratio0 = n0/(s+r0);
        double ratio1 = n1/(s+r1);
        double ratio2 = z2/(s+ 1);
        g = SQR(ratio0) + SQR(ratio1) + SQR(ratio2) - 1;
        if (g > 0) { s0 = s; }
        else if (g < 0) { s1 = s; }
        //else { std::cout << i << std::endl; break; }
        else { break; }
    }
    return s;
}

// For Ellipse
double Ellipsoid::RobustLength(const double v0, const double v1) const {
    const double v0abs = abs(v0);
    const double v1abs = abs(v1);
    return std::max(v0abs, v1abs) * sqrt(1 + SQR(std::min(v0abs, v1abs)/std::max(v0abs, v1abs)));
}

// For Ellipse
double Ellipsoid::GetRoot(const double r0, const double z0, const double z1, double g) const {
    const double n0 = r0*z0;
    double s0 = z1 - 1;
    double s1 = (g < 0 ? 0 : RobustLength(n0, z1) -1);
    double s = 0;
    int maxIterations = 100;
    for (int i=0; i<maxIterations; ++i) {
        s = (s0 + s1)/2;
        if (s == s0 or s == s1) { break; }
        double ratio0 = n0/(s + r0);
        double ratio1 = z1/(s + 1);
        g = SQR(ratio0) + SQR(ratio1) - 1;
        if (g > 0) { s0 = s; }
        else if (g < 0) { s1 = s; }
        else { break; }
    }
    return s;
}

}
