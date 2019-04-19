/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef SCRIPT_INTERFACE_PAIR_CRITERIA_PAIR_CRITERIA_HPP
#define SCRIPT_INTERFACE_PAIR_CRITERIA_PAIR_CRITERIA_HPP

#include "../auto_parameters/AutoParameters.hpp"
#include "core/pair_criteria/pair_criteria.hpp"
#include <string>

namespace ScriptInterface {
namespace PairCriteria {

class PairCriterion : public AutoParameters<PairCriterion> {
public:
  virtual std::shared_ptr<::PairCriteria::PairCriterion>
  pair_criterion() const = 0;
  Variant call_method(std::string const &method,
                      VariantMap const &parameters) override {
    if (method == "decide") {
      return pair_criterion()->decide(boost::get<int>(parameters.at("id1")),
                                      boost::get<int>(parameters.at("id2")));
    }
    throw std::runtime_error("Unknown method called.");
  }
};

class DistanceCriterion : public PairCriterion {
public:
  DistanceCriterion() : m_c(new ::PairCriteria::DistanceCriterion()) {
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

class EnergyCriterion : public PairCriterion {
public:
  EnergyCriterion() : m_c(new ::PairCriteria::EnergyCriterion()) {
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
  std::shared_ptr<::PairCriteria::EnergyCriterion> m_c;
};

class BondCriterion : public PairCriterion {
public:
  BondCriterion() : m_c(new ::PairCriteria::BondCriterion()) {
    add_parameters(
        {{"bond_type",
          [this](Variant const &v) { m_c->set_bond_type(get_value<int>(v)); },
          [this]() { return m_c->get_bond_type(); }}});
  }

  std::shared_ptr<::PairCriteria::PairCriterion>
  pair_criterion() const override {
    return m_c;
  }

private:
  std::shared_ptr<::PairCriteria::BondCriterion> m_c;
};

} /* namespace PairCriteria */
} /* namespace ScriptInterface */

#endif
