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

/**
 * @file
 * @ref ScafacosContext implements the interface of the ScaFaCoS bridge.
 * It is further derived for the coulombic and dipolar versions of ScaFaCoS.
 */

#include "config/config.hpp"

#if defined(SCAFACOS) or defined(SCAFACOS_DIPOLES)

#include "scafacos/ScafacosContextBase.hpp"

#include <cstddef>
#include <string>
#include <tuple>

namespace detail {
std::tuple<Utils::Vector3d const &, Utils::Vector3i, std::size_t>
get_system_params();
}

template <class ScafacosInterface>
struct ScafacosContext : virtual public ScafacosContextBase,
                         public ScafacosInterface {
  using ScafacosContextBase::ScafacosContextBase;
  using ScafacosInterface::ScafacosInterface;

  void update_system_params() override {
    auto params = detail::get_system_params();
    ScafacosInterface::set_runtime_parameters(
        std::get<0>(params).data(), std::get<1>(params).data(),
        static_cast<int>(std::get<2>(params)));
  }

  std::string get_method() const override {
    return ScafacosInterface::get_method();
  }

  std::string get_parameters() const override {
    return ScafacosInterface::get_parameters();
  }
};

#endif // SCAFACOS or SCAFACOS_DIPOLES
