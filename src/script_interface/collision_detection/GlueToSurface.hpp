/*
 * Copyright (C) 2021-2024 The ESPResSo project
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

#include <config/config.hpp>

#ifdef COLLISION_DETECTION
#ifdef VIRTUAL_SITES_RELATIVE

#include "Protocol.hpp"

#include "core/collision_detection/GlueToSurface.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <variant>

namespace ScriptInterface::CollisionDetection {

class GlueToSurface : public Protocol {
  using CoreClass = ::CollisionDetection::GlueToSurface;

public:
  GlueToSurface() {
    add_parameters({
        {"bond_centers", AutoParameter::read_only,
         [this]() {
           return get_bond_variant_by_id(
               std::get<CoreClass>(*m_protocol).bond_centers);
         }},
        {"distance", AutoParameter::read_only,
         [this]() { return std::get<CoreClass>(*m_protocol).distance; }},
        {"bond_vs", AutoParameter::read_only,
         [this]() {
           return get_bond_variant_by_id(
               std::get<CoreClass>(*m_protocol).bond_vs);
         }},
        {"part_type_vs", AutoParameter::read_only,
         [this]() { return std::get<CoreClass>(*m_protocol).part_type_vs; }},
        {"distance_glued_particle_to_vs", AutoParameter::read_only,
         [this]() {
           return std::get<CoreClass>(*m_protocol).dist_glued_part_to_vs;
         }},
        {"part_type_to_be_glued", AutoParameter::read_only,
         [this]() {
           return std::get<CoreClass>(*m_protocol).part_type_to_be_glued;
         }},
        {"part_type_to_attach_vs_to", AutoParameter::read_only,
         [this]() {
           return std::get<CoreClass>(*m_protocol).part_type_to_attach_vs_to;
         }},
        {"part_type_after_glueing", AutoParameter::read_only,
         [this]() {
           return std::get<CoreClass>(*m_protocol).part_type_after_glueing;
         }},
        {"_bond_centers", AutoParameter::read_only,
         [this]() { return std::get<CoreClass>(*m_protocol).bond_centers; }},
        {"_bond_vs", AutoParameter::read_only,
         [this]() { return std::get<CoreClass>(*m_protocol).bond_vs; }},
    });
  }
  std::shared_ptr<::CollisionDetection::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::CollisionDetection::ActiveProtocol> m_protocol;

protected:
  void do_initialize(VariantMap const &params) override {
    auto const bond_center = find_bond_id(params.contains("_bond_centers")
                                              ? params.at("_bond_centers")
                                              : params.at("bond_centers"));
    auto const bond_vs =
        find_bond_id(params.contains("_bond_vs") ? params.at("_bond_vs")
                                                 : params.at("bond_vs"));
    m_protocol = std::make_shared<::CollisionDetection::ActiveProtocol>(
        CoreClass(get_value<double>(params, "distance"), bond_center, bond_vs,
                  get_value<int>(params, "part_type_vs"),
                  get_value<double>(params, "distance_glued_particle_to_vs"),
                  get_value<int>(params, "part_type_to_be_glued"),
                  get_value<int>(params, "part_type_to_attach_vs_to"),
                  get_value<int>(params, "part_type_after_glueing")));
  }
};

} // namespace ScriptInterface::CollisionDetection

#endif // COLLISION_DETECTION
#endif // VIRTUAL_SITES_RELATIVE
