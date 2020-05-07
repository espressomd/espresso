/*
  Copyright (C) 2015-2020 The ESPResSo project

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

#ifndef ESPRESSO_SCRIPTINTERFACE_GENERICDD_GENERIC_DD_HPP
#define ESPRESSO_SCRIPTINTERFACE_GENERICDD_GENERIC_DD_HPP

#include "script_interface/ScriptInterface.hpp"
#include <generic-dd/generic_dd.hpp>
#include <generic-dd/metric.hpp>
#include <string>

extern int this_node;

namespace ScriptInterface {
namespace GenericDD {

class SIGenericDD : public ScriptInterfaceBase {
public:
  SIGenericDD() = default;

  Variant call_method(const std::string &name,
                      const VariantMap &parameters) override {
    if (name == "repart") {
      std::string m = boost::get<std::string>(parameters.at("metric"));

      ::generic_dd::repartition(m);
    } else if (name == "command") {
      std::string m = boost::get<std::string>(parameters.at("cmd"));
      ::generic_dd::command(m);
    }
    return {};
  }
};

} // namespace GenericDD
} // namespace ScriptInterface

#endif
