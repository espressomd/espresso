/*
 * Copyright (C) 2016-2022 The ESPResSo project
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

/* Unit tests for the ScriptInterface::ParticleSlice class. */

#define BOOST_TEST_MODULE ParticleSlice test
#define BOOST_TEST_DYN_LINK
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/unit_test.hpp>

#include "script_interface/Exception.hpp"
#include <script_interface/Variant.hpp>
#include <script_interface/particle_data/ParticleSlice.hpp>
#include <string>
#include <vector>
using namespace ScriptInterface;

bool validate_exception_msg(const Exception &ex) {
  std::string err_msg{ex.what()};
  return err_msg == "Values must be of type vector.";
}

BOOST_AUTO_TEST_CASE(
    non_vector_type_raises_error_in_SetParticleParametersVisitor) {
  Particles::ParticleSlice p_slice{};
  Variant non_vector_variant{1};
  VariantMap params{{"name", std::string("v")},
                    {std::string("values"), non_vector_variant}};
  std::string param_name{"v"};
  Context *context{nullptr};
  std::vector<int> particle_ids{0};
  std::shared_ptr<CellSystem::CellSystem> cell_structure{nullptr};
  std::shared_ptr<Interactions::BondedInteractions> bonded_ias{nullptr};

  BOOST_CHECK_EXCEPTION(std::visit(
                            [&](auto &&vals) {
                              Particles::SetParticleParametersVisitor{}(
                                  particle_ids, param_name, vals, context,
                                  cell_structure, bonded_ias);
                            },
                            params.at("values")),
                        Exception, validate_exception_msg);
}
