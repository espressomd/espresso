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

#include <exception>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <type_traits>

#include "TclCommand.hpp"
#include "TclScriptInterface.hpp"
#include "utils/NumeratedContainer.hpp"
#include "utils/make_bimap.hpp"

namespace ScriptInterface {
namespace Tcl {

/**
 * @brief Tcl style instance management by id.
 *
 * When registered as a command, new objects can be created by "new <name>
 * <params>",
 * the params are parsed by the newly created object. The command returns an id
 * to
 * identify the new object. It can be deleted from Tcl by "delete <id>" and
 * printed
 * by id. Without parameters the command prints all objects.
 */
class TclScriptInterfaceManager : public TclCommand {
public:
  TclScriptInterfaceManager(
      Tcl_Interp *interp,
      const boost::bimap<std::string, std::string> &name_map)
      : TclCommand(interp), m_name_map(name_map) {}

  virtual void parse_from_string(std::list<std::string> &argv) override {
    switch (argv.size()) {
    case 0:
      /* Print everything */
      Tcl_AppendResult(interp(), print_to_string().c_str(), nullptr);
      break;
    case 1: {
      /* Print object by id */
      const int id = parse_id(argv);

      Tcl_AppendResult(interp(), do_print(id).c_str(), nullptr);

      break;
    }
    default: {
      int id;
      /* The first argument has either 'new' to create a new object,
       * 'delete' to remove one or a valid id to set parameters on a
       * existing object.
       */
      if (argv.front() == "new") {
        /* Pop the 'new' */
        argv.pop_front();
        id = do_new(argv);

      } else if (argv.front() == "delete") {
        /* Pop the 'delete' */
        argv.pop_front();

        const int id = parse_id(argv);
        do_delete(id);

        /* Don't parse parameteres for a deleted object. */
        break;
      } else {
        id = parse_id(argv);
      }

      /* Pass remaining arguments to the
       * object for parsing. */

      m_om[id].parse_from_string(argv);

      std::stringstream ss;
      ss << id;
      Tcl_AppendResult(interp(), ss.str().c_str(), 0);
      break;
    }
    }
  }

  virtual std::string print_to_string() const override {
    std::ostringstream ss;

    /* Format as TCL list of lists. */
    for (auto const &it : m_om) {
      ss << "{ " << do_print(it.first) << " } ";
    }

    return ss.str();
  }

private:
  virtual int do_new(std::list<std::string> &argv) {
    std::cout << __PRETTY_FUNCTION__ << ": argv[0] = " << argv.front()
              << std::endl;
    std::cout << __PRETTY_FUNCTION__
              << ": m_name_map.size() = " << m_name_map.size() << std::endl;

    for (auto it = m_name_map.begin(); it != m_name_map.end(); ++it) {
      std::cout << __PRETTY_FUNCTION__ << ": " << it->left << " -> "
                << it->right << std::endl;
    }

    /* Get the internal name for the class to create.
     * will throw if it does not exists. */
    auto const &class_name = m_name_map.left.at(argv.front());

    std::cout << __PRETTY_FUNCTION__ << ": class_name = " << class_name
              << std::endl;

    /* Pop the name */
    argv.pop_front();

    /* Construct the object */
    auto o = TclScriptInterface(class_name, interp());

    /* Add to list and get an index. */
    return m_om.add(o);
  }

  virtual void do_delete(int id) { m_om.remove(id); }

  virtual std::string do_print(int id) const {
    /* Look up the object to print */
    auto const &o = m_om[id];

    /* Get the tcl name */
    auto const &tcl_name = m_name_map.right.at(o.name());

    /* Print the name and the parameters */
    return tcl_name + " " + o.print_to_string();
  }

  int parse_id(std::list<std::string> &argv) const {
    int id;

    std::stringstream ss(argv.front());
    ss >> id;
    if (ss.fail()) {
      throw std::invalid_argument("Expected integer id");
    }
    /* Pop the id */
    argv.pop_front();

    return id;
  }

  Utils::NumeratedContainer<TclScriptInterface> m_om;
  boost::bimap<std::string, std::string> const &m_name_map;
};

} /* namespace Tcl */
} /* namespace ScriptInterface */

#endif
