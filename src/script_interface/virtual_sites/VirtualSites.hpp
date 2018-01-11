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

#ifndef SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_HPP
#define SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_HPP

#include "auto_parameters/AutoParameters.hpp"
#include "config.hpp"
#include "core/virtual_sites.hpp"

namespace ScriptInterface {
namespace VirtualSites {

#ifdef VIRTUAL_SITES
class VirtualSites : public AutoParameters {
public:
  VirtualSites() {
    add_parameters(
        {{"have_velocity",
          [this](const Variant &v) {
            virtual_sites()->set_have_velocity( // NOLINT, bug in clang-tidy-4.0
                get_value<bool>(v));            // NOLINT, bug in clang-tidy-4.0
          },
          [this]() {
            return virtual_sites() // NOLINT, bug in clang-tidy-4.0
                ->have_velocity(); // NOLINT, bug in clang-tidy-4.0
          }}});
  }
  /** Vs implementation we are wrapping */
  virtual std::shared_ptr<::VirtualSites> virtual_sites() = 0;
};

#endif

} /* namespace VirtualSites */
} /* namespace ScriptInterface */
#endif
