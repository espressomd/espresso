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

#ifndef SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_HPP
#define SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_HPP

#include "core/cluster_analysis/Cluster.hpp"
#include "script_interface/ScriptInterface.hpp"

#include <utils/Factory.hpp>

namespace ScriptInterface {
namespace ClusterAnalysis {

class Cluster : public AutoParameters<Cluster> {
public:
  Cluster() = default;
  Variant call_method(std::string const &method,
                      VariantMap const &parameters) override {
    if (method == "particle_ids") {
      return m_cluster->particles;
    }
    if (method == "size") {
      return (int)m_cluster->particles.size();
    }
    if (method == "longest_distance") {
      return m_cluster->longest_distance();
    }
    if (method == "radius_of_gyration") {
      return m_cluster->radius_of_gyration();
    }
    if (method == "fractal_dimension") {
      double mean_sq_residual;
      double df;
      std::tie(df, mean_sq_residual) =
          m_cluster->fractal_dimension(boost::get<double>(parameters.at("dr")));
      return std::vector<double>({df, mean_sq_residual});
    }
    if (method == "center_of_mass") {
      return m_cluster->center_of_mass();
    }
    return false;
  }
  void set_cluster(std::shared_ptr<::ClusterAnalysis::Cluster> &c) {
    m_cluster = c;
  }

private:
  std::shared_ptr<::ClusterAnalysis::Cluster> m_cluster;
};

} /* namespace ClusterAnalysis */
} /* namespace ScriptInterface */

#endif
