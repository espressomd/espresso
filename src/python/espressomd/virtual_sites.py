# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
from .__init__ import has_features
from .script_interface import ScriptInterfaceHelper, script_interface_register

if has_features("VIRTUAL_SITES"):
    @script_interface_register
    class ActiveVirtualSitesHandle(ScriptInterfaceHelper):

        """Handle for the virtual sites implementation active in the core

        This should not be used directly.

        Attributes
        ----------
        implementation :
            instance of a virtual sites implementation

        """
        _so_name = "VirtualSites::ActiveVirtualSitesHandle"

    @script_interface_register
    class VirtualSitesOff(ScriptInterfaceHelper):

        """Virtual sites implementation which does nothing (default)"""
        _so_name = "VirtualSites::VirtualSitesOff"


if has_features("VIRTUAL_SITES_INERTIALESS_TRACERS"):
    @script_interface_register
    class VirtualSitesInertialessTracers(ScriptInterfaceHelper):

        """Virtual sites which are advected with an lb fluid without inertia.
        Forces are on them are transferred to the fluid instantly.

        """
        _so_name = "VirtualSites::VirtualSitesInertialessTracers"


if has_features("VIRTUAL_SITES_RELATIVE"):
    @script_interface_register
    class VirtualSitesRelative(ScriptInterfaceHelper):

        """Virtual sites implementation placing virtual sites relative to other
        particles. See :ref:`Rigid arrangements of particles` for details.

        """
        _so_name = "VirtualSites::VirtualSitesRelative"
