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

#ifndef __RHOMBOID_HPP
#define __RHOMBOID_HPP

#include "Shape.hpp"

namespace Shapes {
  struct Rhomboid : public Shape {
    virtual const std::string name() const { return std::string("Rhomboid"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);
    Parameters all_parameters() const;
    std::map<std::string, Variant> get_parameters();
    void set_parameter(const std::string &name, const Variant &value);
    
    /** corner of the rhomboid */
    Utils::Vector3d pos;
    /** edges adjacent to the corner */
    Utils::Vector3d a;
    Utils::Vector3d b;
    Utils::Vector3d c;
    /** rhomboid direction. (+1 outside -1 inside interaction direction)*/
    double direction;
  };
}

#endif
