/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#include "initialize.hpp"

#include "ActiveVirtualSitesHandle.hpp"
#include "VirtualSitesInertialessTracers.hpp"
#include "VirtualSitesOff.hpp"
#include "VirtualSitesRelative.hpp"
#include "VirtualSitesCenterOfMass.hpp"

namespace ScriptInterface {
namespace VirtualSites {

void initialize(Utils::Factory<ObjectHandle> *om) {
#ifdef VIRTUAL_SITES
  om->register_new<VirtualSitesOff>("VirtualSites::VirtualSitesOff");
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  om->register_new<VirtualSitesInertialessTracers>(
      "VirtualSites::VirtualSitesInertialessTracers");
#endif
#ifdef VIRTUAL_SITES_RELATIVE
  om->register_new<VirtualSitesRelative>("VirtualSites::VirtualSitesRelative");
#endif
#ifdef VIRTUAL_SITES_CENTER_OF_MASS
  om->register_new<VirtualSitesCenterOfMass>("VirtualSites::VirtualSitesCenterOfMass");
#endif
  om->register_new<ActiveVirtualSitesHandle>(
      "VirtualSites::ActiveVirtualSitesHandle");
#endif
}

} // namespace VirtualSites
} // namespace ScriptInterface
