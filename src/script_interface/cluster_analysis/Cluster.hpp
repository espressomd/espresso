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

#pragma once

#include "core/BoxGeometry.hpp"
#include "core/cluster_analysis/Cluster.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/particle_data/ParticleList.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace ClusterAnalysis {

class Cluster : public AutoParameters<Cluster> {
public:
  Cluster() = default;
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "particle_ids") {
      return m_cluster->particles;
    }
    if (method == "particles") {
      return m_particle_list.lock()->call_method(
          "by_ids", {{"id_selection", m_cluster->particles}});
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
      auto const [df, mean_sq_residual] =
          m_cluster->fractal_dimension(get_value<double>(parameters.at("dr")));
      return std::vector<double>({df, mean_sq_residual});
    }
    if (method == "center_of_mass") {
      return m_cluster->center_of_mass();
    }
    return {};
  }
  void set_cluster(std::shared_ptr<::ClusterAnalysis::Cluster> const &c) {
    m_cluster = c;
  }
  void set_particle_list(std::weak_ptr<Particles::ParticleList> const &handle) {
    m_particle_list = handle;
  }

private:
  std::shared_ptr<::ClusterAnalysis::Cluster> m_cluster;
  std::weak_ptr<Particles::ParticleList> m_particle_list;
};

} /* namespace ClusterAnalysis */
} /* namespace ScriptInterface */
