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

#include "TclScriptInterfaceManager.hpp"

namespace ScriptInterface {
namespace Tcl {

class TclConstraintManager : public TclScriptInterfaceManager {
public:
  TclConstraintManager(Tcl_Interp *interp,
                       const boost::bimap<std::string, std::string> &name_map)
      : TclScriptInterfaceManager(interp, name_map) {}

private:
  /* We have to adapt here because constraints with different shapes are
   * not different types that are created by the factory, but the shape
   * is passed in as a parameter */
  int do_new(std::list<std::string> &argv) override {
    /* Get the internal name for the class to create.
     * will throw if it does not exists. */
    auto const &class_name = m_name_map.left.at(argv.front());

    /* Pop the name */
    argv.pop_front();

    /* Construct the constraint */
    auto constraint = TclScriptInterface("Constraints::Constraint", interp());

    std::cout << __PRETTY_FUNCTION__
              << " constraint = " << constraint.script_object().get()
              << std::endl;

    /* Set the given shape */
    constraint.script_object()->set_parameter("shape", class_name);

    /* Add to list and get an index. */
    return m_om.add(constraint);
  }

  virtual std::string do_print(int id) const {
    /* Look up the object to print */
    auto const &o = m_om[id];

    auto const shape_name =
        boost::get<std::string>(o.script_object()->get_parameter("shape"));

    /* Get the tcl name */
    std::string tcl_name = m_name_map.right.at(shape_name);

    /* Print the name and the parameters */
    auto params = o.script_object()->get_parameters();

    std::stringstream ss;
    for (auto const &p : params) {
      /* Filter out shape */
      if (p.first != "shape")
        ss << " " << p.first << " " << p.second;
    }

    return tcl_name + ss.str();
  }
};

} /* namespace ScriptInterface */
} /* Tcl */
