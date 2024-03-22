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

#include "config/config.hpp"

#include "Caliper.hpp"

#include <script_interface/ObjectHandle.hpp>

#include <utils/Factory.hpp>

namespace ScriptInterface::Profiler {

void initialize(Utils::Factory<ObjectHandle> *om) {
#ifdef CALIPER
  om->register_new<Caliper>("ScriptInterface::Profiler::Caliper");
#else
  static_cast<void>(om);
#endif
}

} // namespace ScriptInterface::Profiler
