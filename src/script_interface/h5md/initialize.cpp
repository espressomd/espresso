/*
  Copyright (C) 2015-2018 The ESPResSo project

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

#include "config.hpp"
#ifdef H5MD
#include "h5md.hpp"
#include "initialize.hpp"
#include "io/writer/h5md_core.hpp"

namespace ScriptInterface {
namespace Writer {
void initialize() {
  ScriptInterface::register_new<ScriptInterface::Writer::H5mdScript>(
      "ScriptInterface::Writer::H5mdScript");
}
} /* namespace Writer */
} /* namespace ScriptInterface */
#endif
