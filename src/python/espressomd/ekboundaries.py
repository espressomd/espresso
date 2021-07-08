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
from .script_interface import ScriptObjectRegistry, ScriptInterfaceHelper, script_interface_register
from .__init__ import has_features

# TODO: feature-def

if has_features(["LB_BOUNDARIES"]):
    @script_interface_register
    class EKBoundaries(ScriptObjectRegistry):

        """
        Creates a set of EK boundaries.

        """

        _so_name = "EKBoundaries::EKBoundaries"

        def add(self, *args, **kwargs):
            """
            Adds a boundary to the set.
            Either a valid boundary is an argument,
            or a valid set of parameters to create a boundary.

            """

            if len(args) == 1:
                if isinstance(args[0], EKBoundary):
                    ekboundary = args[0]
                else:
                    raise TypeError(
                        "Either a EKBoundary object or key-value pairs for the parameters of a EKBoundary object need to be passed.")
            else:
                ekboundary = EKBoundary(**kwargs)
            self.call_method("add", object=ekboundary)
            return ekboundary

        def remove(self, ekboundary):
            """
            Removes a boundary from the set.

            Parameters
            ----------
            ekboundary : :obj:`EKBoundary`
                         The boundary to be removed from the set.

            """

            self.call_method("remove", object=ekboundary)

        def clear(self):
            """
            Removes all boundaries.

            """

            self.call_method("clear")

        def size(self):
            return self.call_method("size")

        def empty(self):

            return self.call_method("empty")

    @script_interface_register
    class EKBoundary(ScriptInterfaceHelper):

        """
        Creates a EK boundary.

        """

        _so_name = "EKBoundaries::EKBoundary"
