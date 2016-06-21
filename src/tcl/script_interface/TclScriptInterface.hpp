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

#include <memory>
#include <string>

#include "utils/Factory.hpp"

#include "TclCommand.hpp"
#include "script_interface/ScriptInterface.hpp"

namespace ScriptInterface {
namespace Tcl {

/**
 * @brief Tcl interface to @c ScriptInterfaceBase.
 *
 * This class encapsulates a ScripterfaceBase for access
 * from Tcl. It owns the corresponding C++ object.
 * Copies are shallow, so all copies point to the same
 * C++ object. New objects are only created if the
 * explicit contructor is used.
 */
class TclScriptInterface : public TclCommand {
public:
  typedef typename Utils::Factory<ScriptInterfaceBase> Factory;

  TclScriptInterface(std::string const &name, Tcl_Interp *interp)
      : TclCommand(interp),
        /* Create a new c++ object */
        m_so(Factory::make(name)) {}

  /**
   * @brief Print parameters in Tcl formatting.
   */
  std::string print_to_string() const override;

  /**
   * @brief Parse arguments from list of strings.
   */
  void parse_from_string(std::list<std::string> &argv) override;

  const std::string name() const { return m_so->name(); }

private:
  std::shared_ptr<ScriptInterfaceBase> m_so;
};
}
} /* namespace */

#endif
