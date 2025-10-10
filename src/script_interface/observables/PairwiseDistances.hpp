/*
 * Copyright (C) 2025 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PAIRWISEDISTANCES_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PAIRWISEDISTANCES_HPP

#include "core/observables/PairwiseDistances.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include <memory>
#include <vector>

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class PairwiseDistances
    : public AutoParameters<PairwiseDistances<CoreObs>, Observable> {
  using Base = AutoParameters<PairwiseDistances<CoreObs>, Observable>;

public:
  using Base::Base;

  PairwiseDistances() {
    this->add_parameters({{"ids", AutoParameter::read_only,
                           [this]() { return m_observable->ids(); }},
                          {"target_ids", AutoParameter::read_only,
                           [this]() { return m_observable->target_ids(); }}});
  }

  void do_construct(VariantMap const &params) override {
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_observable =
          make_shared_from_args<CoreObs, std::vector<int>, std::vector<int>>(
              params, "ids", "target_ids");
    });
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} // namespace Observables
} // namespace ScriptInterface

#endif
