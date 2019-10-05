/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __SLITPORE_HPP
#define __SLITPORE_HPP

#include "Shape.hpp"

namespace Shapes {
class Slitpore : public Shape {
public:
  Slitpore()
      : m_pore_mouth(0.0), m_upper_smoothing_radius(0.0),
        m_lower_smoothing_radius(0.0), m_channel_width(0.0), m_pore_width(0.0),
        m_pore_length(0.0), m_dividing_plane(0.0) {}

  void calculate_dist(const Utils::Vector3d &pos, double &dist,
                      Utils::Vector3d &vec) const override;

  double &pore_mouth() { return m_pore_mouth; }
  double &upper_smoothing_radius() { return m_upper_smoothing_radius; }
  double &lower_smoothing_radius() { return m_lower_smoothing_radius; }
  double &channel_width() { return m_channel_width; }
  double &pore_width() { return m_pore_width; }
  double &pore_length() { return m_pore_length; }
  double &dividing_plane() { return m_dividing_plane; }
  const double &dividing_plane() const { return m_dividing_plane; }

private:
  double m_pore_mouth;
  double m_upper_smoothing_radius;
  double m_lower_smoothing_radius;
  double m_channel_width;
  double m_pore_width;
  double m_pore_length;
  double m_dividing_plane;
};
} // namespace Shapes

#endif
