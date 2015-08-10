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

#ifndef __SHAPE_HPP
#define __SHAPE_HPP

#include <string>
#include <iostream>
#include "grid.hpp"

using namespace std;

namespace Shapes {
  struct Shape {
    virtual int calculate_dist(const double *ppos, double *dist, double *vec) = 0;
    /* Human readable name of the shape. */
    static std::string name() { return std::string("Shape"); }
  };

  struct Wall : public Shape {
    static std::string name() { return std::string("Wall"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);  

    /** normal vector on the plane. */
    double n[3];
    /** distance of the wall from the origin. */
    double d;
  };

  struct Sphere : public Shape {
    std::string name() { return std::string("Sphere"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    double pos[3];
    double rad;
    double direction;
  };


  struct Cylinder : public Shape {
    static std::string name() { return std::string("Cylinder"); }
    virtual int calculate_dist(const double *ppos, double *dist, double *vec);

    /** center of the cylinder. */
    double pos[3];
    /** Axis of the cylinder .*/
    double axis[3];
    /** cylinder radius. */
    double rad;
    /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
    double length;
    /** cylinder direction. (+1 outside -1 inside interaction direction)*/
    double direction;
  };

  struct SpheroCylinder : public Cylinder {
    static std::string name() { return std::string("Spherocylinder"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);
  };

  struct Pore : public Shape {
    static std::string name() { return std::string("Pore"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    /** center of the cylinder. */
    double pos[3];
    /** Axis of the cylinder .*/
    double axis[3];
    /** cylinder radius. */
    double rad_left;
    double rad_right;
    double smoothing_radius;
    /** cylinder length. (!!!NOTE this is only the half length of the cylinder.)*/
    double length;
    double outer_rad_left;
    double outer_rad_right;
  };

  struct Slitpore : public Shape {
    static std::string name() { return std::string("Slitpore"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    double pore_mouth;
    double upper_smoothing_radius;
    double lower_smoothing_radius;
    double channel_width;
    double pore_width;
    double pore_length;
  };

  struct Maze : public Shape {
    static std::string name() { return std::string("Maze"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    /** number of spheres. */
    double nsphere;
    /** dimension of the maze. */
    double dim;
    /** sphere radius. */
    double sphrad;
    /** cylinder (connecting the spheres) radius*/
    double cylrad;
  };

  struct Stomatocyte : public Shape {
    static std::string name() { return std::string("Stomatocyte"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    /** Stomatocyte position. */
    double position_x;
    double position_y;
    double position_z;

    /** Stomatocyte orientation. */
    double orientation_x;
    double orientation_y;
    double orientation_z;

    /** Stomatocyte dimensions. */
    double outer_radius;
    double inner_radius;
    double layer_width;

    /** Inside/Outside (+1 outside -1 inside interaction direction)*/
    double direction;
  };

  struct HollowCone : public Shape {
    static std::string name() { return std::string("HollowCone"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    /** Hollow cone position. */
    double position_x;
    double position_y;
    double position_z;

    /** Hollow cone orientation. */
    double orientation_x;
    double orientation_y;
    double orientation_z;

    /** Hollow cone dimensions. */
    double outer_radius;
    double inner_radius;
    double width;
    double opening_angle;

    /** Inside/Outside (+1 outside -1 inside interaction direction)*/
    double direction;
  };

  struct Rhomboid : public Shape {
    static std::string name() { return std::string("Rhomboid"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);
    /** corner of the rhomboid */
    double pos[3];
    /** edges adjacent to the corner */
    double a[3];
    double b[3];
    double c[3];
    /** rhomboid direction. (+1 outside -1 inside interaction direction)*/
    double direction;
  };

  struct Plane : public Shape {
    static std::string name() { return std::string("Plane"); }
    int calculate_dist(const double *ppos, double *dist, double *vec);

    double pos[3];
  };

}

#endif
