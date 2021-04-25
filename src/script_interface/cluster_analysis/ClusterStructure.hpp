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

#ifndef SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP
#define SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP

#include "core/cluster_analysis/ClusterStructure.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/cluster_analysis/Cluster.hpp"
#include "script_interface/pair_criteria/pair_criteria.hpp"

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
            };
          },
          [this]() { return m_pc; }}});
  };
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "get_cluster") {
      // Note: Cluster objects are generated on the fly, to avoid having to
      // store a script interface object for all clusters (which can by
      // thousands)
      auto c = std::dynamic_pointer_cast<Cluster>(
          context()->make_shared("ClusterAnalysis::Cluster", {}));
      c->set_cluster(
          m_cluster_structure.clusters.at(get_value<int>(parameters.at("id"))));

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
};

} /* namespace ClusterAnalysis */
} /* namespace ScriptInterface */

#endif
