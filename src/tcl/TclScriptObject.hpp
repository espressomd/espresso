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

#ifndef __TCLSCRIPTOBJECT_HPP
#define __TCLSCRIPTOBJECT_HPP

/** Tcl interface to ScriptObjects
 */

#include "TclCommand.hpp"
#include "script_interface/ScriptObject.hpp"

class TclScriptObject : public TclCommand {
public:
  TclScriptObject(ScriptInterface::ScriptObject &so, Tcl_Interp* _interp) : m_so(so), TclCommand(_interp) {};
  std::string print_to_string();
  void parse_from_string(std::list<std::string> &argv);
  using TclCommand::create_command;
  void create_command() {
    create_command(m_so.name());
  }
  const std::string name() { return m_so.name(); };

private:
  ScriptInterface::ScriptObject &m_so;
};

#endif
