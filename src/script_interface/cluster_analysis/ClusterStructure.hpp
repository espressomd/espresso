/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP
#define SCRIPT_INTERFACE_CLUSTER_ANALYSIS_CLUSTER_STRUCTURE_HPP

#include "ScriptInterface.hpp"
#include "core/cluster_analysis.hpp"
#include "core/utils/Factory.hpp"
#include "../pair_criteria/pair_criteria.hpp"


namespace ScriptInterface {
namespace ClusterAnalysis {

class ClusterStructure : public AutoParameters {
public:
  ClusterStructure(){
    add_parameters({
                    {"pair_criterion",
                     [this](Variant const &value) {
                       m_pc =
                           get_value<std::shared_ptr<PairCriteria::PairCriterion>>(value);
                       if (m_pc) {
                         m_cluster_structure->set_pair_criterion(m_pc->pair_criterion());
                       };

                     },
                     [this]() {
                       return (m_pc != nullptr) ? m_pc->id() : ObjectId();
                     }}});
  };
  const std::string name() const override { return "ClusterAnalysis::ClusterStructure"; }
  virtual Variant call_method(std::string const &method,
                              VariantMap const &parameters) override {
    if (method == "get_cluster") {
      
      // Note: Cluster objects are generated on the fly, to avoid having to store
      // a script interface object for all clusters (which can by thousands)
      std::shared_ptr<Cluster> c = std::make_shared<Cluster>();
      c->set_cluster(m_cluster_structure->clusters.at(boost::get<int>(parameters.at("id"))));
      return Variant(c->id());
    }
    if (method == "cluster_ids") {
      std::vector<int> cluster_ids;
      for (auto it : m_cluster_structure->clusters) {
         cluster_ids.push_back(it.first);
      }
      return cluster_ids;
    }
    if (method == "cid_for_particle") {
      return m_cluster_structure->cluster_id.at(boost::get<int>(parameters.at("pid")));
    }
    if ("method"=="clear") {
      m_cluster_structure->clear();
    }
    if ("name" == "run_for_all_pairs") {
      m_cluster_structure->run_for_all_pairs();
    }
    if ("name" == "run_for_bonded_particles") {
      m_cluster_structure->run_for_bonded_particles();
    }
  }                              
private:
  std::shared_ptr<::ClusterStructure> m_cluster_structure;
  std::shared_ptr<PairCriteria::PairCriterion> m_pc;
};

} /* namespace ClusterAnalysis */
} /* namespace ScriptInterface */

#endif
