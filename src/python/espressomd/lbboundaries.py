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


if any(has_features(i) for i in ["LB_BOUNDARIES", "LB_BOUNDARIES_GPU"]):
    @script_interface_register
    class LBBoundaries(ScriptObjectRegistry):

        """
        Creates a set of lattice-Boltzmann boundaries.

        """

        _so_name = "LBBoundaries::LBBoundaries"

        def add(self, *args, **kwargs):
            """
            Adds a boundary to the set of boundaries.
            Either pass a valid boundary as argument,
            or a valid set of parameters to create a boundary.

            """

            if len(args) == 1:
                if isinstance(args[0], LBBoundary):
                    lbboundary = args[0]
                else:
                    raise TypeError(
                        "Either a LBBoundary object or key-value pairs for the parameters of a LBBoundary object need to be passed.")
            else:
                lbboundary = LBBoundary(**kwargs)
            self.call_method("add", object=lbboundary)
            return lbboundary

        def remove(self, lbboundary):
            """
            Removes a boundary from the set.

            Parameters
            ----------
            lbboundary : :obj:`LBBoundary`
                The boundary to be removed from the set.

            """

            self.call_method("remove", object=lbboundary)

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
    class LBBoundary(ScriptInterfaceHelper):

        """
        Creates a LB boundary from a shape.

        The fluid velocity is limited to :math:`v_{\\mathrm{max}} = 0.20`
        (see *quasi-incompressible limit* in :cite:`kruger17a`,
        chapter 7, page 272), which corresponds to Mach 0.35.

        The relative error in the fluid density between a compressible fluid
        and an incompressible fluid at Mach 0.30 is less than 5% (see
        *constant density assumption* in :cite:`kundu01a` chapter 16, page
        663). Since the speed of sound is :math:`c_s = 1 / \\sqrt{3}` in LB
        velocity units in a D3Q19 lattice, the velocity limit at Mach 0.30
        is :math:`v_{\\mathrm{max}} = 0.30 / \\sqrt{3} \\approx 0.17`.
        At Mach 0.35 the relative error is around 6% and
        :math:`v_{\\mathrm{max}} = 0.35 / \\sqrt{3} \\approx 0.20`.

        Parameters
        ----------
        shape : :obj:`espressomd.shapes.Shape`
            The shape from which to build the boundary.
        velocity : (3,) array_like of :obj:`float`, optional
            The boundary slip velocity. By default, a velocity of zero is used
            (no-slip boundary condition).

        """

        _so_name = "LBBoundaries::LBBoundary"
        _so_bind_methods = ("get_force",)
