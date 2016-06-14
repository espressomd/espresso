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

#ifndef TCL_TCLSCRIPTOBJECT_HPP
#define TCL_TCLSCRIPTOBJECT_HPP

/** Tcl interface to ScriptObjects
 */

#include "TclCommand.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface { namespace Tcl {

class TclScriptInterface : public TclCommand {
public:
  TclScriptInterface(ScriptInterfaceBase &so, Tcl_Interp* interp)
      : m_so(so),
        TclCommand(interp)
  {};

  /**
   * @brief Print parameters in Tcl formatting.
   */
  std::string print_to_string() const override;
  
  /**
   * @brief Parse arguments from list of strings.
   */
  void parse_from_string(std::list<std::string> &argv) override;

  /**
   * @brief Register the command with Tcl.
   */
  void create_command() {
    TclCommand::create_command(m_so.name());
  }
  
  std::string name() const {
    return m_so.name();
  }

private:
  ScriptInterfaceBase &m_so;
};

}} /* namespace */

#endif
