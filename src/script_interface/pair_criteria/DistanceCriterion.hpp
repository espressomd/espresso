/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_PAIR_CRITERIA_DISTANCE_CRITERION_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_PAIR_CRITERIA_DISTANCE_CRITERION_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/pair_criteria/PairCriterion.hpp"

#include "core/pair_criteria/DistanceCriterion.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace PairCriteria {

class DistanceCriterion : public PairCriterion {
public:
  DistanceCriterion()
      : m_c(std::make_shared<::PairCriteria::DistanceCriterion>()) {
    add_parameters(
        {{"cut_off",
          [this](Variant const &v) { m_c->set_cut_off(get_value<double>(v)); },
          [this]() { return m_c->get_cut_off(); }}});
  }

  std::shared_ptr<::PairCriteria::PairCriterion>
  pair_criterion() const override {
    return m_c;
  }

private:
  std::shared_ptr<::PairCriteria::DistanceCriterion> m_c;
};

} /* namespace PairCriteria */
} /* namespace ScriptInterface */

#endif
