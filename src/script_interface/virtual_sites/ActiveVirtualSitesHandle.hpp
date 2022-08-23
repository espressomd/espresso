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

#ifndef SCRIPT_INTERFACE_VIRTUAL_SITES_ACTIVE_VIRTUAL_SITES_HANDLE_HPP
#define SCRIPT_INTERFACE_VIRTUAL_SITES_ACTIVE_VIRTUAL_SITES_HANDLE_HPP

#include "config/config.hpp"

#ifdef VIRTUAL_SITES

#include "VirtualSites.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "core/errorhandling.hpp"
#include "core/virtual_sites.hpp"

#include <memory>

namespace ScriptInterface {
namespace VirtualSites {

class ActiveVirtualSitesHandle
    : public AutoParameters<ActiveVirtualSitesHandle> {
public:
  ActiveVirtualSitesHandle() {
    add_parameters({{"implementation",
                     [this](Variant const &value) {
                       m_active_implementation =
                           get_value<std::shared_ptr<VirtualSites>>(value);
                       ::set_virtual_sites(
                           m_active_implementation->virtual_sites());
                     },
                     [this]() { return m_active_implementation; }}});
  }

private:
  std::shared_ptr<VirtualSites> m_active_implementation;
};
} /* namespace VirtualSites */
} /* namespace ScriptInterface */
#endif // VIRTUAL_SITES
#endif
