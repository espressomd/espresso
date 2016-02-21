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

#include "TclConstraintManager.hpp"
#include "TclScriptObject.hpp"
#include "TclScriptObjectManager.hpp"
#include "constraints/ConstraintList.hpp"

#include <iostream>

#ifdef CONSTRAINTS

using namespace Utils;

namespace Constraints {
namespace Tcl {

std::string ConstraintManager::print_one(int id) {
  std::ostringstream os;
  os << m_names[id] << " " << TclScriptObject(*m_objects[id], interp).print_to_string();

  return os.str();
}


const int ConstraintManager::get_id(const std::string &s) const {
  std::stringstream ss(s);
  int id;
  ss >> id;
  if(ss.fail()) {
    throw std::string("usage()");
  }

  return id;
}

void ConstraintManager::parse_from_string(std::list<std::string> &argv) {

  if(argv.front() == "delete") {
    argv.pop_front();
    const int id = get_id(argv.front());
    argv.pop_front();
    
    auto c = m_objects[id];    
    Constraints::list.remove_constraint(c);
    
    m_shapes.erase(id);
    m_objects.remove(id);
    m_names.erase(id);
    
    return;
  }

  if(argv.front() == "mindist_position") {
    argv.pop_front();
    argv.pop_front();
    argv.pop_front();

    std::stringstream ss;
    ss << 0;
  
    Tcl_AppendResult(interp, ss.str().c_str(), 0);

    return;
  }
  
  /** print one constraint */
  if(argv.size() == 1) {
    const int id = get_id(argv.front());
    print_one(id);

    return;
  }
  /** The tcl name of the constraint */
  const std::string name = argv.front();
  argv.pop_front();

  int id;
  
  /** Check if this is plain constraint */
  if(ManagedNumeratedContainer<Constraint>::factory_type::has_builder(name)) {
    std::cout << name << std::endl;
    id = m_objects.add(name);
    TclScriptObject(*m_objects[id], interp).parse_from_string(argv);

    /** If it's not one of the other types it is a GeometryConstraint */
  } else {
    std::cout << name << std::endl;
    id = m_objects.add("geometry_constraint");
    auto c = m_objects[id] ;
    auto s = Shapes::Factory::make(name);

    TclScriptObject(*s, interp).parse_from_string(argv);        
    TclScriptObject(*c, interp).parse_from_string(argv);
    
    c->set_parameter("shape", s->get_id());

    m_shapes[id] = s;
  }

  std::cout << m_objects[id]->name() << std::endl;
  Constraints::list.add_constraint(m_objects[id]);
  m_names[id] = name;
  
  std::stringstream ss;
  ss << id;
  
  Tcl_AppendResult(interp, ss.str().c_str(), 0);
}

std::string ConstraintManager::print_to_string() {  
  std::ostringstream ss;

  for(auto &it: m_objects) {
    ss << "{ " << print_one(it.first) << " } ";
  }

  return ss.str();    
}

}
}

#endif /* CONSTRAINTS */
