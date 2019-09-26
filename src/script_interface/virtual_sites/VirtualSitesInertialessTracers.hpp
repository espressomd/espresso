/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP
#define SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP

#include "VirtualSites.hpp"
#include "config.hpp"
#include "core/virtual_sites/VirtualSitesInertialessTracers.hpp"
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

namespace ScriptInterface {
namespace VirtualSites {

class VirtualSitesInertialessTracers : public VirtualSites {
public:
  VirtualSitesInertialessTracers()
      : m_virtual_sites(new ::VirtualSitesInertialessTracers()){};
  /** Vs implementation we are wrapping */
  std::shared_ptr<::VirtualSites> virtual_sites() override {
    return m_virtual_sites;
  };

private:
  std::shared_ptr<::VirtualSitesInertialessTracers> m_virtual_sites;
};

} /* namespace VirtualSites */
} /* namespace ScriptInterface */
#endif
#endif
