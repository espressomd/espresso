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

#ifndef SCRIPT_INTERFACE_RDF_HPP
#define SCRIPT_INTERFACE_RDF_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/observables/RDF.hpp"

#include <boost/range/algorithm.hpp>

#include <cstddef>
#include <iterator>
#include <memory>
#include <vector>

namespace ScriptInterface {
namespace Observables {

class RDF : public AutoParameters<RDF, Observable> {

public:
  RDF() {
    this->add_parameters(
        {{"ids1", AutoParameter::read_only,
          [this]() { return rdf_observable()->ids1(); }},
         {"ids2", AutoParameter::read_only,
          [this]() { return rdf_observable()->ids2(); }},
         {"n_r_bins", AutoParameter::read_only,
          [this]() { return static_cast<int>(rdf_observable()->n_r_bins); }},
         {"min_r", AutoParameter::read_only,
          [this]() { return rdf_observable()->min_r; }},
         {"max_r", AutoParameter::read_only,
          [this]() { return rdf_observable()->max_r; }}});
  }

  void do_construct(VariantMap const &params) override {
    m_observable = make_shared_from_args<::Observables::RDF, std::vector<int>,
                                         std::vector<int>, int, double, double>(
        params, "ids1", "ids2", "n_r_bins", "min_r", "max_r");
  }

  std::shared_ptr<::Observables::RDF> rdf_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<::Observables::RDF> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
