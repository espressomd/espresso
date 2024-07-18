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
#include "core/cluster_analysis/ClusterStructure.hpp"
#include "core/system/System.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/cluster_analysis/Cluster.hpp"
#include "script_interface/pair_criteria/PairCriterion.hpp"
#include "script_interface/particle_data/ParticleList.hpp"
#include "script_interface/system/System.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace ClusterAnalysis {

class ClusterStructure : public AutoParameters<ClusterStructure> {
public:
  ClusterStructure() : m_pc(nullptr) {
    add_parameters(
        {{"pair_criterion",
          [this](Variant const &value) {
            m_pc =
                get_value<std::shared_ptr<PairCriteria::PairCriterion>>(value);
            if (m_pc) {
              m_cluster_structure.set_pair_criterion(m_pc->pair_criterion());
            }
          },
          [this]() { return m_pc; }}});
  }

  void do_construct(VariantMap const &params) override {
    auto local_params = params;
    local_params.erase("system");
    ObjectHandle::do_construct(local_params);
    auto system_si =
        get_value<std::shared_ptr<System::System>>(params, "system");
    m_particle_list = get_value<std::shared_ptr<Particles::ParticleList>>(
        system_si->get_parameter("part"));
    m_cluster_structure.attach(system_si->get_system().box_geo);
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "get_cluster") {
      auto const cluster_id = get_value<int>(parameters.at("id"));
      // Note: Cluster objects are generated on the fly, to avoid having to
      // store a script interface object for all clusters (which can be
      // thousands)
      auto c = std::dynamic_pointer_cast<Cluster>(
          context()->make_shared("ClusterAnalysis::Cluster", {}));
      c->set_cluster(m_cluster_structure.clusters.at(cluster_id));
      c->set_particle_list(m_particle_list);

      return c;
    }
    if (method == "cluster_ids") {
      std::vector<int> cluster_ids;
      for (const auto &it : m_cluster_structure.clusters) {
        cluster_ids.push_back(it.first);
      }
      return cluster_ids;
    }
    if (method == "n_clusters") {
      return int(m_cluster_structure.clusters.size());
    }
    if (method == "cid_for_particle") {
      return m_cluster_structure.cluster_id.at(
          get_value<int>(parameters.at("pid")));
    }
    if (method == "clear") {
      m_cluster_structure.clear();
    }
    if (method == "run_for_all_pairs") {
      m_cluster_structure.run_for_all_pairs();
    }
    if (method == "run_for_bonded_particles") {
      m_cluster_structure.run_for_bonded_particles();
    }
    return {};
  }

private:
  ::ClusterAnalysis::ClusterStructure m_cluster_structure;
  std::shared_ptr<PairCriteria::PairCriterion> m_pc;
  std::weak_ptr<Particles::ParticleList> m_particle_list;
};

} /* namespace ClusterAnalysis */
} /* namespace ScriptInterface */
