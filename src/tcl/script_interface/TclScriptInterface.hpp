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

/*
 * Make operator<< work with the variant. The visitor
 * that is actually invoked by boost lives in the boost
 * namespace and finds these functions. The other possibility
 * would be to define them in std and rely on ADL, but this is not legal
 * because std::vector<T> is not a user defined type.
 * (Works, though)
 */

#warning move out of namespace std;
namespace std {

/* Source: http://stackoverflow.com/a/6693088/3198615 */
template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  std::ostringstream ss;

  std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(ss, " "));
  ss << v.back();

  out << ss.str();

  return out;
}

template <typename T, int n>
std::ostream &operator<<(std::ostream &out, const Vector<n, T> &v) {
  std::ostringstream ss;

  std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(ss, " "));
  ss << v.back();

  out << ss.str();

  return out;
}
}

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
        m_so(ScriptInterfaceBase::make_shared(name)) {}

  /**
   * @brief Print parameters in Tcl formatting.
   */
  std::string print_to_string() const override;

  /**
   * @brief Parse arguments from list of strings.
   */
  void parse_from_string(std::list<std::string> &argv) override;

  const std::string name() const { return m_so->name(); }
  std::shared_ptr<ScriptInterfaceBase> script_object() { return m_so; }
  const std::shared_ptr<ScriptInterfaceBase> script_object() const {
    return m_so;
  }

private:
  std::shared_ptr<ScriptInterfaceBase> m_so;
};
}
} /* namespace */

#endif
