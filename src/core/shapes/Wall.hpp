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

#ifndef __WALL_HPP
#define __WALL_HPP

#include "Shape.hpp"
#include "Vector.hpp"

namespace Shapes {
  struct Wall : public Shape {
    virtual const std::string name() const { return std::string("Wall"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);  
    Parameters &all_parameters() const;
    Parameters get_parameters();
    void set_parameter(const std::string &name, const Variant &value);
    
    /** normal vector on the plane. */
    Vector3d n;
    /** distance of the wall from the origin. */
    double d;
  };
}

#endif
