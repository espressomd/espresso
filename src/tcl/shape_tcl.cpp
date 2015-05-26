/*
  Copyright (C) 2010,2011,2012,2013,2014,2015 The ESPResSo project
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

#include "parser.hpp"
#include "Shape.hpp"
#include "shape_tcl.hpp"
#include <list>
#include <string>

int tclcommand_shape_parse_wall(Tcl_Interp *interp, std::list<std::string> &args, Shape::Wall *s) {
  s->n[0] = s->n[1] = s->n[2] = 0.0;
  s->d = 0.0;
  while(!args.empty()) {
    if(args.front() == "normal") {
      args.pop_front();
      for(int i = 0; i < 3; i++) {
        if(!STR_IS_D(args.front(), s->n[i])) {
          Tcl_AppendResult(interp, "constraint wall normal <nx> <ny> <nz> expected", (char *) NULL);
          return (TCL_ERROR);
        } else {
          args.pop_front();
        }
      }
    }
    else if(args.front() == "dist") {
      args.pop_front();
      if(!STR_IS_D(args.front(), s->d)) {
        Tcl_AppendResult(interp, "constraint wall dist <d> expected", (char *) NULL);
        return (TCL_ERROR);
      } else
        args.pop_front();
    }
    else
      break;
  }
  double norm = SQR(s->n[0])+SQR(s->n[1])+SQR(s->n[2]);
  if (norm < 1e-10) {
    Tcl_AppendResult(interp, "usage: wall normal <nx> <ny> <nz> dist <d>",
		     (char *) NULL);
    return (TCL_ERROR);    
  }

  return TCL_OK;

}

int tclcommand_shape_parse_sphere(Tcl_Interp *interp, std::list<std::string> &args, Shape::Sphere *s) {
  s->pos[0] = 0.0;
  s->pos[1] = 0.0;
  s->pos[2] = 0.0;
  s->rad = 0;
  s->direction = -1;

  while(!args.empty()) {
    if(args.front() == "center") {
      args.pop_front();
      for(int i = 0; i < 3; i++) {
        if(args.empty() || !STR_IS_D(args.front(), s->pos[i])) {
          Tcl_AppendResult(interp, "constraint sphere center <x> <y> <z> expected", (char *) NULL);
          return (TCL_ERROR);
        } else {
          args.pop_front();
        }
      }
    } else if(args.front() == "radius") {
      args.pop_front();
      for(int i = 0; i < 3; i++) {
        if(args.empty() || !STR_IS_D(args.front(), s->rad)) {
          Tcl_AppendResult(interp, "constraint sphere radius <r> expected", (char *) NULL);
          return (TCL_ERROR);
        } else {
          args.pop_front();
        }
      }
    } else if(args.front() == "direction") {
      args.pop_front();
      if (args.front() == "inside")
        s->direction = -1;
      else if (args.front() == "outside")
        s->direction = 1;
      else if (Tcl_GetDouble(interp, args.front().c_str(), &(s->direction)) == TCL_ERROR) {
        Tcl_AppendResult(interp, "-1/1 or inside/outside is expected", (char *) NULL);
        return (TCL_ERROR);
      }
      args.pop_front();
    }
  }

  return TCL_OK;
}

int tclcommand_shape_parse_cylinder(Tcl_Interp *interp, std::list<std::string> &args, Shape::Cylinder *s) {
  s->pos[0] = 0.0;
  s->pos[1] = 0.0;
  s->pos[2] = 0.0;
  s->axis[0] = 0.0;
  s->axis[1] = 0.0;
  s->axis[2] = 0.0;
  s->rad = 0;
  s->length = 0;
  s->direction = 0;
 
  while(!args.empty()) {
    if(args.front() == "center") {
      args.pop_front();
      for(int i = 0; i < 3; i++) {
        if(args.empty() || !STR_IS_D(args.front(), s->pos[i])) {
          Tcl_AppendResult(interp, "constraint sphere center <x> <y> <z> expected", (char *) NULL);
          return (TCL_ERROR);
        } else {
          args.pop_front();
        }
      }
    }
    return TCL_OK;
  }
}
int tclcommand_shape_parse_rhomboid(Tcl_Interp *interp, std::list<std::string> &args, Shape::Rhomboid *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_shape_parse_maze(Tcl_Interp *interp, std::list<std::string> &args, Shape::Maze *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_shape_parse_pore(Tcl_Interp *interp, std::list<std::string> &args, Shape::Pore *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_shape_parse_slitpore(Tcl_Interp *interp, std::list<std::string> &args, Shape::Slitpore *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_shape_parse_stomatocyte(Tcl_Interp *interp, std::list<std::string> &args, Shape::Stomatocyte *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_shape_parse_hollow_cone(Tcl_Interp *interp, std::list<std::string> &args, Shape::HollowCone *s) {
  while(!args.empty()) {
    ;
  }
  return TCL_OK;
}

int tclcommand_parse_shape(Tcl_Interp *interp, std::list<std::string> &args, Shape::Shape **s) {
  std::string name = args.front();
  args.pop_front();

  if(name == "wall") {
    *s = new Shape::Wall;
    return tclcommand_shape_parse_wall(interp, args, dynamic_cast<Shape::Wall *>(*s));
  }
  else if(name == "sphere") {
    *s = new Shape::Sphere;
    return tclcommand_shape_parse_sphere(interp, args, dynamic_cast<Shape::Sphere *>(*s));
  }
  else if(name == "cylinder") {
    *s = new Shape::Cylinder;
    return tclcommand_shape_parse_cylinder(interp, args, dynamic_cast<Shape::Cylinder *>(*s));
  }
  else if(name == "rhomboid") {
    *s = new Shape::Rhomboid;
    return tclcommand_shape_parse_rhomboid(interp, args, dynamic_cast<Shape::Rhomboid *>(*s));
  }
  else if(name == "maze") {
    *s = new Shape::Maze;
    return tclcommand_shape_parse_maze(interp, args, dynamic_cast<Shape::Maze *>(*s));
  }
  else if(name == "pore") {
    *s = new Shape::Pore;
    return tclcommand_shape_parse_pore(interp, args, dynamic_cast<Shape::Pore *>(*s));
  }
  else if(name == "slitpore") {
    *s = new Shape::Slitpore;
    return tclcommand_shape_parse_slitpore(interp, args, dynamic_cast<Shape::Slitpore *>(*s));
  }
  else if(name == "stomatocyte") {
    *s = new Shape::Stomatocyte;
    return tclcommand_shape_parse_stomatocyte(interp, args, dynamic_cast<Shape::Stomatocyte *>(*s));
  }
  else if(name == "hollow_cone") {
    *s = new Shape::HollowCone;
    return tclcommand_shape_parse_hollow_cone(interp, args, dynamic_cast<Shape::HollowCone *>(*s));
  }
  else {
    Tcl_AppendResult(interp, "unknown shape '", name.c_str(), "'.", (char *) NULL);
    return TCL_OK;
  }
}
