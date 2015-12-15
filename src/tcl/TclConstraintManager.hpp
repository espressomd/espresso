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
#ifndef CONSTRAINT_TCL_H
#define CONSTRAINT_TCL_H

#include "config.hpp"

#ifdef CONSTRAINTS

#include "parser.hpp"
#include "TclCommand.hpp"
#include "constraints/Constraint.hpp"
#include "ObjectContainer.hpp"

namespace Constraints {
namespace Tcl {

class ConstraintManager : public TclCommand {
 public:
  ConstraintManager(Tcl_Interp *interp) : TclCommand(interp) {}
  void parse_from_string(std::list<std::string> &argv);
  std::string print_to_string();
 private:
  ObjectContainer<Constraints::Constraint> m_objects;
  std::map<int, std::string> m_names;
  std::string print_one(int id);
  const int get_id(const std::string &s) const;
};

}
}

#endif /* CONSTRAINTS */

#endif /* CONSTRAINT_TCL_H */
