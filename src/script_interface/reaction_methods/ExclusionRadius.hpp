/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "core/reaction_methods/ExclusionRadius.hpp"

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace ReactionMethods {

class ExclusionRadius : public AutoParameters<ExclusionRadius> {
  std::shared_ptr<::ExclusionRadius> m_obj;

public:
  ExclusionRadius() {
    add_parameters({{"search_algorithm",
                     [this](Variant const &v) {
                       context()->parallel_try_catch([&]() {
                         auto const key = get_value<std::string>(v);
                         if (key == "order_n") {
                           m_obj->neighbor_search_order_n = true;
                         } else if (key == "parallel") {
                           m_obj->neighbor_search_order_n = false;
                         } else {
                           throw std::invalid_argument(
                               "Unknown search algorithm '" + key + "'");
                         }
                       });
                     },
                     [this]() {
                       if (m_obj->neighbor_search_order_n) {
                         return std::string("order_n");
                       }
                       return std::string("parallel");
                     }},
                    {"exclusion_range",
                     [this](Variant const &v) {
                       context()->parallel_try_catch([&]() {
                         m_obj->set_exclusion_range(get_value<double>(v));
                       });
                     },
                     [this]() { return m_obj->exclusion_range; }},
                    {"exclusion_radius_per_type",
                     [this](Variant const &v) {
                       context()->parallel_try_catch([&]() {
                         m_obj->set_exclusion_radius_per_type(
                             get_value<::ExclusionRadius::map_type>(v));
                       });
                     },
                     [this]() {
                       return make_unordered_map_of_variants(
                           m_obj->exclusion_radius_per_type);
                     }}});
  }

  void do_construct(VariantMap const &params) override {
    context()->parallel_try_catch([&]() {
      m_obj = std::make_shared<::ExclusionRadius>(context()->get_comm());
    });
    if (params.count("exclusion_range")) {
      do_set_parameter("exclusion_range", params.at("exclusion_range"));
    }
    if (params.count("exclusion_radius_per_type")) {
      do_set_parameter("exclusion_radius_per_type",
                       params.at("exclusion_radius_per_type"));
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "check_exclusion_range") {
      auto const pid = get_value<int>(params, "pid");
      if (params.count("ptype")) {
        auto const ptype = get_value<int>(params, "ptype");
        return m_obj->check_exclusion_range(pid, ptype);
      }
      return m_obj->check_exclusion_range(pid);
    }
    return {};
  }

  auto get_instance() const { return m_obj; }
};

} // namespace ReactionMethods
} // namespace ScriptInterface
