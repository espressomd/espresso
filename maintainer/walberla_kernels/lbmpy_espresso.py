#
# Copyright (C) 2021-2023 The ESPResSo project
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
#

import pystencils as ps

import lbmpy.boundaries
import lbmpy.custom_code_nodes

import lbmpy_walberla.additional_data_handler


class BounceBackSlipVelocityUBB(
        lbmpy_walberla.additional_data_handler.UBBAdditionalDataHandler):
    '''
    Dynamic UBB that implements the bounce-back method with slip velocity.
    '''

    def data_initialisation(self, direction):
        '''
        Modified ``indexVector`` initialiser. The "classical" dynamic UBB
        uses the velocity callback as a velocity flow profile generator.
        Here we use that callback as a bounce-back slip velocity generator.
        This way, the dynamic UBB can be used to implement a LB boundary.
        '''
        code = super().data_initialisation(direction)
        assert "InitialisationAdditionalData" in code
        assert "elementInitialiser" in code
        assert "element.vel_0" in code
        bb_vec = self.stencil_info[direction][1]
        cell_args = [f"it.{direction}() + {bb_vec[i]}".replace('+ -', '-')
                     for i, direction in enumerate("xyz")]
        code = [
            "auto const InitialisationAdditionalData = elementInitialiser(",
            f"Cell({', '.join(cell_args)}), blocks, *block);",
            "element.vel_0 = InitialisationAdditionalData[0];",
            "element.vel_1 = InitialisationAdditionalData[1];",
            "element.vel_2 = InitialisationAdditionalData[2];",
        ]
        return "\n".join(code)


class UBB(lbmpy.boundaries.UBB):
    '''
    Velocity bounce back boundary condition, enforcing specified velocity at
    obstacle. This is a patched version of ``lbmpy.boundaries.UBB``, which
    currently doesn't support the bounce back scheme we need.
    '''

    def __call__(self, f_out, f_in, dir_symbol,
                 inv_dir, lb_method, index_field):
        '''
        Modify the assignments such that the source and target pdfs are swapped.
        '''
        assignments = super().__call__(
            f_out, f_in, dir_symbol, inv_dir, lb_method, index_field)

        assert len(assignments) > 0

        out = []
        if len(assignments) > 1:
            out.extend(assignments[:-1])

        neighbor_offset = lbmpy.custom_code_nodes.NeighbourOffsetArrays.neighbour_offset(
            dir_symbol, lb_method.stencil)

        assignment = assignments[-1]
        assert assignment.lhs.field == f_in
        out.append(ps.Assignment(assignment.lhs.get_shifted(*neighbor_offset),
                                 assignment.rhs - f_out(dir_symbol) + f_in(dir_symbol)))
        return out
