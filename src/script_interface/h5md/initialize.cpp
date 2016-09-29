/*
  Copyright (C) 2015,2016 The ESPResSo project

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

#include "initialize.hpp"
#include "ParallelScriptInterface.hpp"
#include "h5md_core.hpp"
#include "h5md.hpp"


namespace ScriptInterface {
namespace Writer {
namespace H5md {
    void initialize() {
        ParallelScriptInterface<ScriptInterface::Writer::H5md::H5md>::register_new(
            "Writer::H5md::File");
    }
} /* namespace H5md */
} /* namespace Writer */
} /* namespace ScriptInterface */

