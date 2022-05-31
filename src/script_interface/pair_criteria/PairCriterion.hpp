/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_PAIR_CRITERIA_PAIR_CRITERION_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_PAIR_CRITERIA_PAIR_CRITERION_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"

#include "core/pair_criteria/PairCriterion.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace PairCriteria {

class PairCriterion : public AutoParameters<PairCriterion> {
public:
  virtual std::shared_ptr<::PairCriteria::PairCriterion>
  pair_criterion() const = 0;
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "decide") {
      return pair_criterion()->decide(get_value<int>(parameters.at("id1")),
                                      get_value<int>(parameters.at("id2")));
    }
    throw std::runtime_error("Unknown method called.");
  }
};

} /* namespace PairCriteria */
} /* namespace ScriptInterface */

#endif
