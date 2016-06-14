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
#ifndef TCL_TCLSCRIPTOBJECTMANAGER_HPP
#define TCL_TCLSCRIPTOBJECTMANAGER_HPP

#include <type_traits>
#include <iostream>
#include <string>
#include <exception>
#include <memory>
#include <list>

#include "TclCommand.hpp"
#include "TclScriptInterface.hpp"
#include "utils/NumeratedContainer.hpp"
#include "utils/Factory.hpp"

namespace ScriptInterface { namespace Tcl {

/** @brief Tcl style instance management by id.
 *
 * When registered as a command, new objects can be created by "new <name> <params>",
 * the params are parsed by the newly created object. The command returns an id to
 * identify the new object. It can be deleted from Tcl by "delete <id>" and printed
 * by id. Without parameters the command prints all objects.
 */
template<class T, class Factory = Utils::Factory<T> >
class TclScriptInterfaceManager : public TclCommand {
 public:
  TclScriptInterfaceManager(Tcl_Interp *_interp, std::map<std::string, std::string> const& names) : TclCommand(_interp) {
    static_assert(std::is_base_of<ScriptInterfaceBase, T>::value, "Type has to be subclass of ScriptInterfaceBase");
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

            auto o = Factory::make(argv.front());
            
            id = m_om.add(std::move(o));
            argv.pop_front();
          } else {
            std::stringstream ss(argv.front());
            ss >> id;
            if(ss.fail()) {
              throw std::invalid_argument("Expected integer id");
            }
            argv.pop_front();
          }
        
          TclScriptInterfaceBase(*m_om[id], interp).parse_from_string(argv);
          std::stringstream ss;
          ss << id;
          Tcl_AppendResult(interp, ss.str().c_str(), 0);
          break;
        }
    }      
  }

  std::string print_to_string() {
    std::ostringstream ss;
    
    for(auto & it: m_om) {
      ss << "{ " << print_one(it.first) << " } ";
    }

    return ss.str();    
  }


 private:
  Utils::NumeratedContainer<std::unique_ptr<TclScriptInterface> > m_om;

  std::string print_one(const int id) {
    return std::string(m_om.name(id)).append(" ").append(TclScriptInterface(*m_om[id], interp).print_to_string());
  }
};

} /* namespace Tcl */ } /* namespace ScriptInterface */

#endif
