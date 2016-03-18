/*
  Copyright (C) 2015,2016 The ESPResSo project
  
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
#ifndef __TCLSCRIPTOBJECTMANAGER_HPP
#define __TCLSCRIPTOBJECTMANAGER_HPP

#include <type_traits>
#include <iostream>
#include <string>
#include <exception>
#include <memory>
#include <list>

#include "TclCommand.hpp"
#include "script_interface/ScriptObject.hpp"
#include "utils/ManagedNumeratedContainer.hpp"

namespace Tcl {

/** Tcl Interface for a ObjectManager<ScriptObject> */

template<class T>
class TclScriptObjectManager : public TclCommand {
public:
  TclScriptObjectManager(Tcl_Interp *_interp) : TclCommand(_interp) {
    static_assert(std::is_base_of<ScriptInterface::ScriptObject, T>::value, "Type has to be subclass of ScriptObject");
  }

  void parse_from_string(std::list<std::string> &argv) {
    switch(argv.size()) {
    case 0:
      Tcl_AppendResult(interp, print_to_string().c_str(), 0);
      break;
    case 1:
      {
        int id;
        std::stringstream ss(argv.front());
        
        ss >> id;
        if(ss.fail()) {
          throw std::invalid_argument("Expected integer id");
        }        
        argv.pop_front();
        
        Tcl_AppendResult(interp, print_one(id).c_str(), 0);
        
        break;
      }
    default:
      {
        int id;
        if(argv.front() == "new") {
          argv.pop_front();
          id = m_om.add(argv.front());
          argv.pop_front();
        } else {
          std::stringstream ss(argv.front());
          ss >> id;
          if(ss.fail()) {
            throw std::invalid_argument("Expected integer id");
          }
          argv.pop_front();
        }
        
        TclScriptObject(*m_om[id], interp).parse_from_string(argv);
        std::stringstream ss;
        ss << id;
        Tcl_AppendResult(interp, ss.str().c_str(), 0);
        break;
      }
    }      
  }

  std::string print_to_string() {
    std::ostringstream ss;
    
    for(typename Utils::ManagedNumeratedContainer<T>::iterator it = m_om.begin(); it != m_om.end(); ++it) {
      ss << "{ " << print_one(it->first) << " } ";
    }

    return ss.str();    
  }


private:
  Utils::ManagedNumeratedContainer<T> m_om;

  std::string print_one(const int id) {
    return std::string(m_om.name(id)).append(" ").append(TclScriptObject(*m_om[id], interp).print_to_string());
  }
};

}

#endif
