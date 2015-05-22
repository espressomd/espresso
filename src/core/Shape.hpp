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

struct Shape {
  virtual int calculate_dist(const double *ppos, double *dist, double *vec) = 0;
  /* Human readable name of the shape. */
  static std::string name() { return std::string("Shape::"); }
};

struct Wall : public Shape {
  Wall(double *_n, double _d) : d(_d) {
    n[0] = _n[0];
    n[1] = _n[1];
    n[2] = _n[2];
  }
  static std::string name() { return std::string("Wall"); }
  int calculate_dist(const double *ppos, double *dist, double *vec);  

  /** normal vector on the plane. */
  double n[3];
  /** distance of the wall from the origin. */
  double d;
};

struct Sphere : public Shape {
  Sphere(double *_pos, double _rad, double _direction) :
    rad(_rad), direction(_direction) {
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
  }
  std::string name() { return std::string("Sphere"); }
  int calculate_dist(const double *ppos, double *dist, double *vec);

  double pos[3];
  double rad;
  double direction;
};


struct Cylinder : public Shape {
  Cylinder(double *_pos, double *_axis, double _rad, double _length, double _direction) :
    rad(_rad), length(_length), direction(_direction) {
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];
  }
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
  SpheroCylinder(double *_pos, double *_axis, double _rad, double _length, double _direction) : Cylinder(_pos, _axis, _rad, _length, _direction) {}
  static std::string name() { return std::string("spherocylinder"); }
  int calculate_dist(const double *ppos, double *dist, double *vec);
};

struct Pore : public Shape {
  Pore(double *_pos, double *_axis, double _rad_left, double _rad_right, double _smoothing_radius, double _length, double _direction, double _outer_rad_left, double _outer_rad_right) :
    rad_left(_rad_left), rad_right(_rad_right), smoothing_radius(_smoothing_radius), length(_length), outer_rad_left(_outer_rad_left), outer_rad_right(_outer_rad_right)
  {
    pos[0] = _pos[0];
    pos[1] = _pos[1];
    pos[2] = _pos[2];
    axis[0] = _axis[0];
    axis[1] = _axis[1];
    axis[2] = _axis[2];
  }
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
  Slitpore(double _pore_mouth,
           double _upper_smoothing_radius,
           double _lower_smoothing_radius,
           double _channel_width,
           double _pore_width,
           double _pore_length
           ) : pore_mouth(_pore_mouth), upper_smoothing_radius(_upper_smoothing_radius), lower_smoothing_radius(_lower_smoothing_radius), channel_width(_channel_width), pore_width(_pore_width), pore_length(_pore_length) {};
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
  Maze(double _nsphere, double _dim, double _sphrad, double _cylrad) : nsphere(_nsphere), dim(_dim), sphrad(_sphrad), cylrad(_cylrad) {}
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
  Stomatocyte(double _position_x,
              double _position_y,
              double _position_z,
              double _orientation_x,
              double _orientation_y,
              double _orientation_z,
              double _outer_radius,
              double _inner_radius,
              double _layer_width,
              double _direction) :
    position_x(_position_x), position_y(_position_y), position_z(_position_z),
    orientation_x(_orientation_x), orientation_y(_orientation_y), orientation_z(_orientation_z),
    outer_radius(_outer_radius), inner_radius(_inner_radius), layer_width(_layer_width), direction(_direction) {}
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
  HollowCone(double _position_x,
             double _position_y,
             double _position_z,
             double _orientation_x,
             double _orientation_y,
             double _orientation_z,
             double _outer_radius,
             double _inner_radius,
             double _width,
             double _opening_angle,
             double _direction) :
    position_x(_position_x), position_y(_position_y), position_z(_position_z),
    orientation_x(_orientation_x), orientation_y(_orientation_y), orientation_z(_orientation_z),
    outer_radius(_outer_radius), inner_radius(_inner_radius), width(_width), opening_angle(_opening_angle), direction(_direction) {}
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

struct Box : public Shape {
  Box(double _value) : value(_value) {}
  static std::string name() { return std::string("Box"); }
  int calculate_dist(const double *ppos, double *dist, double *vec);

  double value;
};

struct Rhomboid : public Shape {
  Rhomboid(double *_pos, double *_a, double *_b, double *_c, double _direction) : direction(_direction) {
    pos[0] =_pos[0];
    pos[1] =_pos[1];
    pos[2] =_pos[2];
    a[0] = _a[0];
    a[1] = _a[1];
    a[2] = _a[2];
    b[0] = _b[0];
    b[1] = _b[1];
    b[2] = _b[2];
    c[0] = _c[0];
    c[1] = _c[1];
    c[2] = _c[2];
  }
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
  Plane(double *_pos) {
    pos[0] =_pos[0];
    pos[1] =_pos[1];
    pos[2] =_pos[2];
  }
  static std::string name() { return std::string("Plane"); }
  int calculate_dist(const double *ppos, double *dist, double *vec);

  double pos[3];
};

#endif
