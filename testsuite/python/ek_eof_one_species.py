# Copyright (C) 2011-2019 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import pathlib

import sys
import math
import numpy as np
try:
    import vtk
    from vtk.util import numpy_support as VN
    skipIfMissingPythonPackage = utx.no_skip
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module vtk not available, skipping test!")


import espressomd
import espressomd.electrokinetics
import espressomd.shapes
import ek_common

##########################################################################
#                          Set up the System                             #
##########################################################################
# Set the slit pore geometry. The width is the non-periodic part of the
# geometry. The padding is used to ensure that there is no field outside
# the slit.

params_base = dict([
    ('dt', 1.0 / 7),
    ('integration_length', 2300),
    ('agrid', 1. / 3),
    ('density_water', 26.15),
    ('friction', 1.9),
    ('width', 20.0),
    ('thickness', 3.0),
    ('sigma', -0.04),
    ('padding', 6.0),
    ('force', 0.07),
    ('temperature', 1.1),
    ('viscosity_kinematic', 1.7),
    ('bjerrum_length', 0.8),
    ('sigma', -0.04),
    ('valency', 1.0),
])
params_base['density_counterions'] = -2.0 * \
    params_base['sigma'] / params_base['width']

axis = "@TEST_SUFFIX@"
params = {
    "x": dict([
        ('box_x', params_base['thickness']),
        ('box_y', params_base['thickness']),
        ('box_z', params_base['width'] + 2 * params_base['padding']),
        ('ext_force_density', [params_base['force'], 0.0, 0.0]),
        ('wall_normal_1', [0, 0, 1]),
        ('wall_normal_2', [0, 0, -1]),
        ('periodic_dirs', (0, 1)),
        ('non_periodic_dir', 2),
        ('n_roll_index', 0),
        ('calculated_pressure_xy', 0.0),
        ('calculated_pressure_yz', 0.0)
    ]),
    "y": dict([
        ('box_x', params_base['width'] + 2 * params_base['padding']),
        ('box_y', params_base['thickness']),
        ('box_z', params_base['thickness']),
        ('ext_force_density', [0.0, params_base['force'], 0.0]),
        ('wall_normal_1', [1, 0, 0]),
        ('wall_normal_2', [-1, 0, 0]),
        ('periodic_dirs', (1, 2)),
        ('non_periodic_dir', 0),
        ('n_roll_index', 1),
        ('calculated_pressure_xz', 0.0),
        ('calculated_pressure_yz', 0.0)
    ]),
    "z": dict([
        ('box_x', params_base['thickness']),
        ('box_y', params_base['width'] + 2 * params_base['padding']),
        ('box_z', params_base['thickness']),
        ('ext_force_density', [0.0, 0.0, params_base['force']]),
        ('wall_normal_1', [0, 1, 0]),
        ('wall_normal_2', [0, -1, 0]),
        ('periodic_dirs', (0, 2)),
        ('non_periodic_dir', 1),
        ('n_roll_index', 2),
        ('calculated_pressure_xy', 0.0),
        ('calculated_pressure_xz', 0.0)
    ])
}[axis]


def bisection():
    # initial parameters for bisection scheme
    size = math.pi / (2.0 * params_base['width'])
    pnt0 = 0.0
    pntm = pnt0 + size
    pnt1 = pnt0 + 1.9 * size

    # the bisection scheme
    tol = 1.0e-08
    while size > tol:
        val0 = ek_common.solve(
            pnt0,
            params_base['width'],
            params_base['bjerrum_length'],
            params_base['sigma'],
            params_base['valency'])
        val1 = ek_common.solve(
            pnt1,
            params_base['width'],
            params_base['bjerrum_length'],
            params_base['sigma'],
            params_base['valency'])
        valm = ek_common.solve(
            pntm,
            params_base['width'],
            params_base['bjerrum_length'],
            params_base['sigma'],
            params_base['valency'])

        if (val0 < 0.0 and val1 > 0.0):
            if valm < 0.0:
                pnt0 = pntm
                size = size / 2.0
                pntm = pnt0 + size
            else:
                pnt1 = pntm
                size = size / 2.0
                pntm = pnt1 - size
        elif (val0 > 0.0 and val1 < 0.0):
            if valm < 0.0:
                pnt1 = pntm
                size = size / 2.0
                pntm = pnt1 - size
            else:
                pnt0 = pntm
                size = size / 2.0
                pntm = pnt0 + size
        else:
            sys.exit("Bisection method fails:\n"
                     "Tuning of domain boundaries may be required.")
    return pntm


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["ELECTROKINETICS", "EK_BOUNDARIES"])
class ek_eof_one_species(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    xi = bisection()

    def parse_vtk(self, filepath, name, shape):
        reader = vtk.vtkStructuredPointsReader()
        reader.SetFileName(filepath)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()
        points = data.GetPointData()

        return VN.vtk_to_numpy(points.GetArray(name)).reshape(shape, order='F')

    @classmethod
    def setUpClass(cls):
        system = cls.system
        system.box_l = [params['box_x'], params['box_y'], params['box_z']]
        system.time_step = params_base['dt']
        system.cell_system.skin = 0.1
        system.thermostat.turn_off()

        # Set up the (LB) electrokinetics fluid
        ek = cls.ek = espressomd.electrokinetics.Electrokinetics(
            agrid=params_base['agrid'],
            lb_density=params_base['density_water'],
            viscosity=params_base['viscosity_kinematic'],
            friction=params_base['friction'],
            T=params_base['temperature'],
            prefactor=params_base['bjerrum_length'] *
            params_base['temperature'],
            stencil="linkcentered")

        counterions = cls.counterions = espressomd.electrokinetics.Species(
            density=params_base['density_counterions'],
            D=0.3,
            valency=params_base['valency'],
            ext_force_density=params['ext_force_density'])
        ek.add_species(counterions)

        # Set up the walls confining the fluid and carrying charge
        ek_wall1 = espressomd.ekboundaries.EKBoundary(
            charge_density=params_base['sigma'] /
            params_base['padding'],
            shape=espressomd.shapes.Wall(
                normal=params['wall_normal_1'],
                dist=params_base['padding']))
        ek_wall2 = espressomd.ekboundaries.EKBoundary(
            charge_density=params_base['sigma'] /
            params_base['padding'],
            shape=espressomd.shapes.Wall(
                normal=params['wall_normal_2'],
                dist=-(params_base['padding'] + params_base['width'])))
        system.ekboundaries.add(ek_wall1)
        system.ekboundaries.add(ek_wall2)
        system.actors.add(ek)

        # Integrate the system
        system.integrator.run(params_base['integration_length'])

    def test(self):
        # compare the various quantities to the analytic results
        total_velocity_difference = 0.0
        total_density_difference = 0.0
        total_pressure_difference_xx = 0.0
        total_pressure_difference_yy = 0.0
        total_pressure_difference_zz = 0.0
        total_pressure_difference_xy = 0.0
        total_pressure_difference_yz = 0.0
        total_pressure_difference_xz = 0.0

        system = self.system
        ek = self.ek
        counterions = self.counterions
        for i in range(
                int(system.box_l[params['non_periodic_dir']] / params_base['agrid'])):
            if (i *
                params_base['agrid'] >= params_base['padding'] and i *
                params_base['agrid'] < system.box_l[params['non_periodic_dir']] -
                    params_base['padding']):
                position = i * params_base['agrid'] - params_base['padding'] - \
                    params_base['width'] / 2.0 + params_base['agrid'] / 2.0

                # density
                index = np.array([int(system.box_l[params['periodic_dirs'][0]] /
                                      (2 * params_base['agrid'])),
                                  int(system.box_l[params['periodic_dirs'][1]] /
                                      (2 * params_base['agrid'])), i])
                index = np.roll(index, params['n_roll_index'])
                measured_density = counterions[index].density
                calculated_density = ek_common.density(
                    position, self.xi, params_base['bjerrum_length'])
                density_difference = abs(measured_density - calculated_density)
                total_density_difference += density_difference

                # velocity
                measured_velocity = ek[index].velocity[int(
                    np.nonzero(params['ext_force_density'])[0])]
                calculated_velocity = ek_common.velocity(
                    position,
                    self.xi,
                    params_base['width'],
                    params_base['bjerrum_length'],
                    params_base['force'],
                    params_base['viscosity_kinematic'],
                    params_base['density_water'])
                velocity_difference = abs(
                    measured_velocity - calculated_velocity)
                total_velocity_difference = total_velocity_difference + \
                    velocity_difference

                # diagonal pressure tensor
                measured_pressure_xx = ek[index].pressure_tensor[(0, 0)]
                calculated_pressure_xx = ek_common.hydrostatic_pressure(
                    ek,
                    (0, 0),
                    system.box_l[params['periodic_dirs'][0]],
                    system.box_l[params['periodic_dirs'][1]],
                    params['box_z'],
                    params_base['agrid'])
                measured_pressure_yy = ek[index].pressure_tensor[(1, 1)]
                calculated_pressure_yy = ek_common.hydrostatic_pressure(
                    ek,
                    (1, 1),
                    system.box_l[params['periodic_dirs'][0]],
                    system.box_l[params['periodic_dirs'][1]],
                    params['box_z'],
                    params_base['agrid'])
                measured_pressure_zz = ek[index].pressure_tensor[(2, 2)]
                calculated_pressure_zz = ek_common.hydrostatic_pressure(
                    ek,
                    (2, 2),
                    system.box_l[params['periodic_dirs'][0]],
                    system.box_l[params['periodic_dirs'][1]],
                    params['box_z'],
                    params_base['agrid'])

                pressure_difference_xx = abs(
                    measured_pressure_xx - calculated_pressure_xx)
                pressure_difference_yy = abs(
                    measured_pressure_yy - calculated_pressure_yy)
                pressure_difference_zz = abs(
                    measured_pressure_zz - calculated_pressure_zz)

                total_pressure_difference_xx = total_pressure_difference_xx + \
                    pressure_difference_xx
                total_pressure_difference_yy = total_pressure_difference_yy + \
                    pressure_difference_yy
                total_pressure_difference_zz = total_pressure_difference_zz + \
                    pressure_difference_zz

                calculated_pressure_offdiagonal = ek_common.pressure_tensor_offdiagonal(
                    position, self.xi, params_base['bjerrum_length'], params_base['force'])
                # xy component pressure tensor
                measured_pressure_xy = ek[index].pressure_tensor[(0, 1)]
                calculated_pressure_xy = 0.0
                if 'calculated_pressure_xy' not in params:
                    calculated_pressure_xy = calculated_pressure_offdiagonal
                pressure_difference_xy = abs(
                    measured_pressure_xy - calculated_pressure_xy)
                total_pressure_difference_xy = total_pressure_difference_xy + \
                    pressure_difference_xy

                # yz component pressure tensor
                measured_pressure_yz = ek[index].pressure_tensor[(1, 2)]
                calculated_pressure_yz = 0.0
                if 'calculated_pressure_yz' not in params:
                    calculated_pressure_yz = calculated_pressure_offdiagonal
                pressure_difference_yz = abs(
                    measured_pressure_yz - calculated_pressure_yz)
                total_pressure_difference_yz = total_pressure_difference_yz + \
                    pressure_difference_yz

                # xz component pressure tensor
                measured_pressure_xz = ek[index].pressure_tensor[(0, 2)]
                calculated_pressure_xz = 0.0
                if 'calculated_pressure_xz' not in params:
                    calculated_pressure_xz = calculated_pressure_offdiagonal
                pressure_difference_xz = abs(
                    measured_pressure_xz - calculated_pressure_xz)
                total_pressure_difference_xz = total_pressure_difference_xz + \
                    pressure_difference_xz

        scale_factor = params_base['agrid'] / params_base['width']
        total_density_difference *= scale_factor
        total_velocity_difference *= scale_factor
        total_pressure_difference_xx *= scale_factor
        total_pressure_difference_yy *= scale_factor
        total_pressure_difference_zz *= scale_factor
        total_pressure_difference_xy *= scale_factor
        total_pressure_difference_yz *= scale_factor
        total_pressure_difference_xz *= scale_factor

        self.assertLess(total_density_difference, 1.0e-04,
                        "Density accuracy not achieved")
        self.assertLess(total_velocity_difference, 1.0e-04,
                        "Velocity accuracy not achieved")
        self.assertLess(total_pressure_difference_xx, 1.0e-04,
                        "Pressure accuracy xx component not achieved")
        self.assertLess(total_pressure_difference_yy, 1.0e-04,
                        "Pressure accuracy yy component not achieved")
        self.assertLess(total_pressure_difference_zz, 1.0e-04,
                        "Pressure accuracy zz component not achieved")
        self.assertLess(total_pressure_difference_xy, 1.0e-04,
                        "Pressure accuracy xy component not achieved")
        self.assertLess(total_pressure_difference_yz, 1.0e-04,
                        "Pressure accuracy yz component not achieved")
        self.assertLess(total_pressure_difference_xz, 1.0e-04,
                        "Pressure accuracy xz component not achieved")

    @skipIfMissingPythonPackage
    def test_vtk(self):
        ek = self.ek
        counterions = self.counterions
        grid_dims = list(
            map(int, np.round(self.system.box_l / params_base['agrid'])))

        # write VTK files
        vtk_root = f"vtk_out/ek_eof_{axis}"
        pathlib.Path(vtk_root).mkdir(parents=True, exist_ok=True)
        path_vtk_boundary = f"{vtk_root}/boundary.vtk"
        path_vtk_velocity = f"{vtk_root}/velocity.vtk"
        path_vtk_potential = f"{vtk_root}/potential.vtk"
        path_vtk_lbdensity = f"{vtk_root}/density.vtk"
        path_vtk_lbforce = f"{vtk_root}/lbforce.vtk"
        path_vtk_density = f"{vtk_root}/lbdensity.vtk"
        path_vtk_flux = f"{vtk_root}/flux.vtk"
        path_vtk_flux_link = f"{vtk_root}/flux_link.vtk"
        if espressomd.has_features('EK_DEBUG'):
            path_vtk_flux_fluc = f"{vtk_root}/flux_fluc.vtk"
        ek.write_vtk_boundary(path_vtk_boundary)
        ek.write_vtk_velocity(path_vtk_velocity)
        ek.write_vtk_potential(path_vtk_potential)
        ek.write_vtk_density(path_vtk_lbdensity)
        ek.write_vtk_lbforce(path_vtk_lbforce)
        counterions.write_vtk_density(path_vtk_density)
        counterions.write_vtk_flux(path_vtk_flux)
        if espressomd.has_features('EK_DEBUG'):
            counterions.write_vtk_flux_fluc(path_vtk_flux_fluc)
        counterions.write_vtk_flux_link(path_vtk_flux_link)

        # load VTK files to check they are correctly formatted
        get_vtk = self.parse_vtk
        vtk_boundary = get_vtk(path_vtk_boundary, "boundary", grid_dims)
        vtk_velocity = get_vtk(path_vtk_velocity, "velocity", grid_dims + [3])
        vtk_potential = get_vtk(path_vtk_potential, "potential", grid_dims)
        vtk_lbdensity = get_vtk(path_vtk_lbdensity, "density_lb", grid_dims)
        get_vtk(path_vtk_lbforce, "lbforce", grid_dims + [3])
        vtk_density = get_vtk(path_vtk_density, "density_1", grid_dims)
        vtk_flux = get_vtk(path_vtk_flux, "flux_1", grid_dims + [3])
        if espressomd.has_features('EK_DEBUG'):
            get_vtk(path_vtk_flux_fluc, "flux_fluc_1", grid_dims + [4])
        get_vtk(path_vtk_flux_link, "flux_link_1", grid_dims + [13])

        # check VTK files against the EK grid
        species_density = np.zeros(grid_dims)
        species_flux = np.zeros(grid_dims + [3])
        ek_potential = np.zeros(grid_dims)
        ek_velocity = np.zeros(grid_dims + [3])
        for i in range(grid_dims[0]):
            for j in range(grid_dims[1]):
                for k in range(grid_dims[2]):
                    index = np.array([i, j, k])
                    species_density[i, j, k] = counterions[index].density
                    species_flux[i, j, k] = counterions[index].flux
                    ek_potential[i, j, k] = ek[index].potential
                    ek_velocity[i, j, k] = ek[index].velocity

        np.testing.assert_allclose(vtk_velocity, ek_velocity, atol=1e-6)
        np.testing.assert_allclose(vtk_potential, ek_potential, atol=1e-6)
        np.testing.assert_allclose(vtk_density, species_density, atol=1e-6)
        np.testing.assert_allclose(vtk_flux, species_flux, atol=1e-6)

        # check VTK files against the EK parameters
        dens = params_base['density_water']
        left_dist = int(params_base['padding'] / params_base['agrid'])
        right_dist = int(-params_base['padding'] / params_base['agrid'])
        thickness = int(params_base['thickness'] / params_base['agrid'])
        i = np.roll([0, 0, right_dist], params['n_roll_index'])
        j = np.roll([thickness, thickness, left_dist], params['n_roll_index'])
        mask_left = np.zeros(grid_dims, dtype=bool)
        mask_left[:j[0], :j[1], :j[2]] = True
        mask_right = np.zeros(grid_dims, dtype=bool)
        mask_right[i[0]:, i[1]:, i[2]:] = True
        mask_outside = np.logical_or(mask_left, mask_right)
        mask_inside = np.logical_not(mask_outside)
        np.testing.assert_allclose(vtk_lbdensity[mask_inside], dens, atol=1e-4)
        np.testing.assert_allclose(vtk_lbdensity[mask_outside], 0, atol=1e-6)
        np.testing.assert_allclose(vtk_boundary[mask_left], 1, atol=1e-6)
        np.testing.assert_allclose(vtk_boundary[mask_left], 1, atol=1e-6)
        np.testing.assert_allclose(vtk_boundary[mask_right], 2, atol=1e-6)
        np.testing.assert_allclose(vtk_boundary[mask_inside], 0, atol=1e-6)


if __name__ == "__main__":
    ut.main()
