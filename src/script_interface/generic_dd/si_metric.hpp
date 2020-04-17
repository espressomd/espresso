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

#ifndef ESPRESSO_SCRIPTINTERFACE_GENERICDD_METRIC_HPP
#define ESPRESSO_SCRIPTINTERFACE_GENERICDD_METRIC_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include <generic-dd/metric.hpp>
#include <string>

namespace ScriptInterface {
namespace GenericDD {

class SIMetric : public AutoParameters<SIMetric> {
public:
  SIMetric() {
    add_parameters({{"metric",
                     [this](const Variant &v) {
                       m_metric_desc = get_value<std::string>(v);
                       m_metric.set_metric(m_metric_desc);
                     },
                     [this]() { return m_metric_desc; }}});
  }

  Variant call_method(const std::string &name,
                      const VariantMap &parameters) override {
    if (name == "average")
      return m_metric.paverage();
    else if (name == "maximum")
      return m_metric.pmax();
    else if (name == "imbalance")
      return m_metric.pimbalance();
    else
      return {};
  }

  const ::generic_dd::repart::Metric &get_metric() const { return m_metric; }

private:
  std::string m_metric_desc;
  ::generic_dd::repart::Metric m_metric;
};

} // namespace GenericDD
} // namespace ScriptInterface

#endif
