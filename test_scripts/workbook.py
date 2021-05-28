# Change working directory from the workspace root to the ipynb file location. Turn this addition
# off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import sys as sysos
import os
import random
import math
import logging
import numpy as np
import argparse
import espressomd
from espressomd.interactions import FeneBond
from espressomd.interactions import AngleHarmonic
from espressomd.io.writer import vtf
from espressomd.magnetostatics import DipolarDirectSumCpu
from espressomd.virtual_sites import VirtualSitesRelative
from espressomd import checkpointing
from espressomd import shapes
from collections.abc import Iterable
from itertools import product
import operator
import numbers
import gzip
import shutil
from datetime import datetime
try:
    import cPickle as pickle
except ImportError:
    import pickle
try:
    sysos.path.append('/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts')
    logging.info(os.getcwd())
except:
    pass
root = logging.getLogger()
root.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sysos.stdout)
handler.setLevel(logging.DEBUG)
root.addHandler(handler)


def args():
    pass


args.alpha = 0.6
args.len = 10
args.h_field = 0
args.eta = 4.1
args.init = (30, 30, 0)
args.slip_v = 1
# %%


def init_filament_params(filament):
    logging.info('init_filament_params runned')
    param_dict = dict(
        num_particles=None, config=[],
        init_pos=(None, None, None),
        type_nonmag=None, type_mag=None, type_virt=None, type_to_be_magnetized=None, fene_k=None,
        fene_r_max=None, fene_r0=None, bond_len=None, ang_seeds=(None, None),
        dip_moments=(None, None, None),
        dip_magnitude=None, magnetize=None, sinTheta=None, cosTheta=None, cosPhi=None, sinPhi=None,
        last_index_used=None,
        realz_indices=None, corner_indices=None, edge_indices=[])
    filament.__dict__.update(param_dict)
    return filament


def set_filament_params(filament):
    logging.info('system_check runned')
    assert core is not None, 'run core() to initialise var necessary ro run'
    logging.info('set_filament_params runned')
    filament.__dict__['type_nonmag'] = 0
    filament.__dict__['type_mag'] = 1
    filament.__dict__['type_virt'] = 2
    filament.__dict__['type_to_be_magnetized'] = 3
    # CAUTION!!!!!!!!!!!!!!!!!!!!!!!
    # type in config doesnt affect the type of real particle in gen_positions(), only affects the length of the chain!
    filament.__dict__['config'] = [0 for x in range(args.len)]
    filament.__dict__['num_particles'] = len(filament.__dict__['config'])
    filament.__dict__['ang_seeds'] = random.random(), random.random()
    # filament.__dict__['sinTheta'] = -1.0 + 2 * filament.__dict__['ang_seeds'][0]
    # filament.__dict__['cosTheta'] = pow((1-pow(filament.__dict__['sinTheta'], 2.)), 0.5)
    # filament.__dict__['cosPhi'] = np.cos(2 * np.pi * filament.__dict__['ang_seeds'][1])
    # filament.__dict__['sinPhi'] = np.sin(2 * np.pi * filament.__dict__['ang_seeds'][1])
    filament.__dict__['sinTheta'] = 0
    filament.__dict__['cosTheta'] = pow((1-pow(filament.__dict__['sinTheta'], 2.)), 0.5)
    bljo = -0.25
    filament.__dict__['cosPhi'] = np.cos(2 * np.pi * bljo)
    filament.__dict__['sinPhi'] = np.sin(2 * np.pi * bljo)
    filament.generator_pos = gen_positions()
    filament.__dict__['bond_len'] = core.sigma+core.sigma*args.alpha
    filament.__dict__['fene_k'] = 40
    filament.__dict__['fene_r0'] = core.sigma*args.alpha
    filament.__dict__['fene_r_max'] = filament.fene_r0*3.
    filament.__dict__['init_pos'] = np.array([x/2 for x in core.box_l])+np.array(args.init)
    print(args.init)
    filament.__dict__['magnetize'] = False
    filament.__dict__['dip_magnitude'] = 9
    if filament.__dict__['magnetize'] == True:
        filament.__dict__['dip_moments'] = 0, 0, 0
    else:
        filament.dip_moments = (
            filament.dip_magnitude * filament.cosTheta * filament.cosPhi, filament.dip_magnitude * filament.cosTheta *
            filament.sinPhi, filament.dip_magnitude * filament.sinTheta)
    filament.__dict__['last_index_used'] = 0
    filament.__dict__['realz_indices'] = []
    filament.__dict__['corner_indices'] = []
    filament.__dict__['edge_indices'] = []
    return filament


def gen_positions():
    def unit_vector(vector):
        return vector / np.linalg.norm(vector)

    def angle_between(v1, v2):
        v1_u, v2_u = unit_vector(v1), unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def generator_real(bond_length, no_of_realz):
        x_move, y_move, z_move = bond_length*cosTheta * \
            cosPhi, bond_length*cosTheta*sinPhi, bond_length*sinTheta
        real_pos_cnt = []
        for n in range(no_of_realz):
            ikser, vajer, zeder = coord_X + n*x_move, coord_Y + n*y_move, coord_Z + n*z_move
            real_pos_cnt.append((ikser, vajer, zeder))
        return real_pos_cnt
    instance = next(globals()[x] for x in globals()
                    if isinstance(globals()[x], Simulation))
    logging.info(sysos.path[-1])
    sinTheta = instance.filament.attributes.sinTheta
    cosTheta = instance.filament.attributes.cosTheta
    cosPhi = instance.filament.attributes.cosPhi
    sinPhi = instance.filament.attributes.sinPhi
    coord_X, coord_Y, coord_Z = instance.filament.attributes.init_pos
    bond_length = instance.filament.attributes.bond_len
    x_move, y_move, z_move = cosTheta * cosPhi, cosTheta*sinPhi, sinTheta
    z, x = np.array((x_move, y_move, z_move)), np.array([0.5, 0.5, 0.5])
    x -= x.dot(z) * z / np.linalg.norm(z)**2
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    range_x, range_y, range_z = [0, ], [0, ], [0, ]
    with open(os.path.join(sysos.path[-1], 'cube_mesh_sigma3_nonula.txt'), mode='r') as source:
        summary = source.readlines()
        for line in summary:
            iks, ipsilon, zed = line.strip().split(sep=' ')
            range_x.append(float(iks))
            range_y.append(float(ipsilon))
            range_z.append(float(zed))
    pass_the_hot_potato = generator_real(
        bond_length, instance.filament.attributes.num_particles)
    grid_of_positions = []
    for (i, j, k) in zip(range_x, range_y, range_z):
        grid_of_positions.append(tuple(map(operator.add, tuple(map(
            operator.add, tuple(np.array(x)*i), tuple(np.array(y)*j))), tuple(np.array(z)*k))))

    for (index, (px, py, pz)) in enumerate(pass_the_hot_potato):
        grid_of_positions_sp = grid_of_positions.copy()
        for (i, (x, y, z)) in enumerate(grid_of_positions_sp):
            grid_of_positions_sp[i] = x+px, y+py, z+pz
        shape_again = np.shape(grid_of_positions_sp)
        yield grid_of_positions_sp, index


def system_check(sys_setter):
    logging.info('system_check runned')
    assert core is not None, 'run core() to initialise var necessary ro run'
    return sys_setter


def init_sys_params(core):
    logging.info('init_sys_params runned')
    param_dict = dict(periodicity=(None, None, None),
                      time_step=None, box_l=(None, None, None),
                      skin=None, sigma=None, virt_sigma=None, min_global_cut=None, temperature=None, gamma=None, seed=None,
                      annealing_steps=None, corner_sigma=None, omega=None)
    core.__dict__.update(param_dict)
    return core


def set_sys_params(core):
    logging.info('set_sys_params runned')
    core.__dict__['periodicity'] = (False, True, True)
    core.__dict__['time_step'] = 0.01
    core.__dict__['box_l'] = (140, 140, 280)
    core.__dict__['sigma'] = 2.0
    core.__dict__['virt_sigma'] = 0.5
    core.__dict__['min_global_cut'] = 0.6
    core.__dict__['skin'] = core.min_global_cut
    core.__dict__['temperature'] = 1
    core.__dict__['gamma'] = 1
    random.seed(os.getpid())
    core.__dict__['seed'] = random.choice(range(os.getpid()))
    core.__dict__['annealing_steps'] = 100
    return core


@set_sys_params
@init_sys_params
def core():
    return core


traceMe = False


def trace(*args):
    if traceMe:
        logging.info('[' + ' '.join(map(str, args)) + ']')


def Private(*privates):
    def onDecorator(aClass):
        class onInstance:
            def __init__(self, *args, **kargs):
                self.wrapped = aClass(*args, **kargs)

            def __getattr__(self, attr):
                trace('get:', attr)  # Others assumed in wrapped
                if attr in privates:
                    raise TypeError('private attribute')
                else:
                    return getattr(self.wrapped, attr)

            def __setattr__(self, attr, value):
                trace('set:', attr, value)
                if attr == 'wrapped':
                    self.__dict__[attr] = value
                elif attr in privates:
                    raise TypeError('private attribute')
                else:
                    setattr(self.wrapped, attr, value)
        return onInstance
    return onDecorator

# %%


@Private('calculate_generic_dipole_zfield_at', 'effective_dipole_field_at_point', 'reset_dipoles')
class Filament():

    def __init__(self):
        logging.info('filament init')

        @set_filament_params
        @init_filament_params
        def filament():
            pass
        self.attributes = filament

    def generator_ids(self, grid_of_points):
        if isinstance(grid_of_points, Iterable):
            shape_again = np.shape(grid_of_points)
            grid_of_indices = np.empty(shape_again[0], dtype=int)
            for (index, element) in enumerate(grid_of_points):
                if index in [1, 5, 21, 25, 100, 104, 120, 124]:
                    grid_of_indices[index] = self.attributes.last_index_used
                    self.attributes.corner_indices.append(
                        self.attributes.last_index_used)
                    self.attributes.last_index_used += 1
                elif index in [11, 51, 70, 110, 15, 55, 74, 114]:
                    grid_of_indices[index] = self.attributes.last_index_used
                    self.attributes.edge_indices.append(
                        self.attributes.last_index_used)
                    self.attributes.last_index_used += 1
                else:
                    grid_of_indices[index] = self.attributes.last_index_used
                    self.attributes.last_index_used += 1
            self.attributes.realz_indices.append(grid_of_indices[0])
            return grid_of_indices

    def generator_types(self, grid_of_indices):
        if isinstance(grid_of_indices, Iterable):
            shape_again = int(np.shape(grid_of_indices)[0])
            grid_of_types = np.empty(shape_again, dtype=int)
            for (index, element) in enumerate(grid_of_indices):
                if index in self.attributes.realz_indices:
                    if self.attributes.magnetize == True:
                        grid_of_types[index] = self.attributes.type_nonmag
                    else:
                        grid_of_types[index] = self.attributes.type_nonmag
                else:
                    grid_of_types[index] = self.attributes.type_virt
            return grid_of_types

    def set_filament(self):
        assert sys is not None, 'You havent initialised the sys!!!'
        self.part_pos = []
        while True:
            try:
                pos, at_real_no = next(self.attributes.generator_pos)
                identity = self.generator_ids(pos)
                ty = self.generator_types(identity)
                self.part_pos.append((identity, ty, pos))
                for (pos, ide, ty) in zip(pos, identity, ty):
                    if ty == 0:
                        sys.part.add(pos=pos, id=ide, type=ty, rotation=(1, 1, 1))
                    if ty == 1:
                        sys.part.add(pos=pos, id=ide, type=ty,
                                     dip=self.attributes.dip_moments, rotation=(1, 1, 1))
                    if ty == 2:
                        sys.part.add(pos=pos, id=ide, type=ty,
                                     rotation=(0, 0, 0))

                for identity in identity:
                    if identity != self.attributes.realz_indices[-1]:
                        espressomd.particle_data.ParticleList(
                        )[identity].vs_auto_relate_to(self.attributes.realz_indices[-1])

            except StopIteration:
                logging.info('last index')
                logging.info(self.attributes.last_index_used)
                logging.info('Particles set!!!')
                break

    def set_rot_intertia_from_units(self):
        for x in self.attributes.realz_indices:
            sys.part[x].rinertia = np.ones(3) * 9.436215089936574
            sys.part[x].mass = 14.154322634904863

    def add_dipolez_to_virts(self):
        komad = [x for x in sys.part.select(
            lambda p: p.type == self.attributes.type_nonmag)]
        for real in komad:
            self.attributes.last_index_used += 1
            sys.part.add(
                pos=real.pos, id=self.attributes.last_index_used, type=self.attributes.type_to_be_magnetized,
                dip=self.attributes.dip_moments,
                rotation=(0, 0, 0))
            espressomd.particle_data.ParticleList(
            )[self.attributes.last_index_used].vs_auto_relate_to(real.id)

    def set_steric(self):
        self.attributes.__dict__.update(dict(steric_inter=True))
        lj_eps = 0.1
        # sys.non_bonded_inter[self.attributes.type_nonmag,
        #                      self.attributes.type_nonmag].wca.set_params(epsilon=lj_eps, sigma=core.sigma)
        # sys.non_bonded_inter[self.attributes.type_mag, self.attributes.type_mag].wca.set_params(
        #     epsilon=lj_eps, sigma=core.sigma)
        sys.non_bonded_inter[self.attributes.type_virt, self.attributes.type_virt].wca.set_params(
            epsilon=lj_eps, sigma=core.virt_sigma)

    def print_sys(self):
        if self.attributes.magnetize == True:
            list_of_stuff = [x for x in sys.part.select(
                lambda p: p.type == self.attributes.type_to_be_magnetized)]
        else:
            list_of_stuff = [x for x in sys.part.select(lambda p: p.type == self.attributes.type_nonmag)]
        for x in list_of_stuff:
            yield x.id, x.type, x.pos, x.quat

    def set_bending_interaction(self):
        angle_harmonic = AngleHarmonic(bend=3.2, phi0=np.pi)
        sys.bonded_inter.add(angle_harmonic)
        list_of_stuff = [x.id for x in sys.part.select(
            lambda p: p.id in self.attributes.realz_indices)]
        for x in list_of_stuff[1:-1]:
            sys.part[x].add_bond((angle_harmonic, x-1, x+1))

    def calculate_generic_dipole_zfield_at(self, dipx, dipy, dipz, mx, my, mz, x, y, z):
        dipole_modulus = self.attributes.dip_magnitude
        if dipole_modulus > 0.:
            dx, dy, dz = x-dipx, y-dipy, z-dipz
            dr = np.sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2))
            dr3, dr5 = pow(dr, 3.0), pow(dr, 5.0)
            mr = dx*mx+dy*my+dz*mz
            Hx, Hy, Hz = 3.*dx*mr/dr5-mx/dr3, 3.*dy*mr/dr5-my/dr3, 3.*dz*mr/dr5-mz/dr3
            return Hx, Hy, Hz
        else:
            return list(0, 0, 0)

    def effective_dipole_field_at_point(self, particle_index, slice):
        storage = [0, 0, 0]
        particle_x, particle_y, particle_z = sys.part[particle_index].pos
        for part_id, part_pos, part_dip in slice:
            if part_id != sys.part[particle_index].id:
                dip_x, dip_y, dip_z = part_pos
                m_dip_x, m_dip_y, m_dip_z = part_dip
                dipole_field = self.calculate_generic_dipole_zfield_at(
                    dip_x, dip_y, dip_z, m_dip_x, m_dip_y, m_dip_z, particle_x, particle_y, particle_z)
                storage[0] = storage[0] + dipole_field[0]
                storage[1] = storage[1] + dipole_field[1]
                storage[2] = storage[2] + dipole_field[2]
        return storage

    def reset_dipoles(self):
        dungeon_witch = sys.part.select(lambda p: p.type == self.attributes.type_to_be_magnetized)
        slice_persists = [(x.id, x.pos, x.dip) for x in dungeon_witch]
        dipole_system_use_to_reset_values = (self.effective_dipole_field_at_point(
            part.id, slice_persists) for part in dungeon_witch)
        H = next(x.H for x in sys.constraints)
        for part in dungeon_witch:
            one = next(dipole_system_use_to_reset_values)
            nula, jedan, dva = H[0]+one[0], H[1]+one[1], H[2]+one[2]
            tri = np.linalg.norm((nula, jedan, dva))
            dipole_x = self.attributes.dip_magnitude*(1.0/np.tanh(
                self.attributes.dip_magnitude*tri)-1.0/(self.attributes.dip_magnitude*tri))*nula/tri
            dipole_y = self.attributes.dip_magnitude*(1.0/np.tanh(
                self.attributes.dip_magnitude*tri)-1.0/(self.attributes.dip_magnitude*tri))*jedan/tri
            dipole_z = self.attributes.dip_magnitude*(1.0/np.tanh(
                self.attributes.dip_magnitude*tri)-1.0/(self.attributes.dip_magnitude*tri))*dva/tri
            part.dip = (dipole_x, dipole_y, dipole_z)

    def magnetize(self):
        assert self.attributes.magnetize, "not in self.attributes.__dict__"
        self.reset_dipoles()

    def excommunicado(self):
        logging.info(self.attributes.realz_indices)
        for element in self.attributes.realz_indices:
            # exclude virt-virt interactions
            slice_master = sys.part.select(
                lambda p: p.vs_relative[0] == element)
            slice_virt = [x.id for x in slice_master if x.type == self.attributes.type_virt]
            for index in slice_virt:
                index_other = slice_virt.copy()
                index_other.remove(index)
                sys.part[index].exclusions = index_other
        # sys.integrator.run(steps=0)
        # print(sys.analysis.energy())

    def bond_cornerz(self):
        self.attributes.bond = espressomd.interactions.FeneBond(
            k=self.attributes.fene_k, d_r_max=self.attributes.fene_r_max, r_0=self.attributes.fene_r0)
        sys.bonded_inter.add(self.attributes.bond)
        for x in range(len(self.attributes.realz_indices[:-1])):
            slice_master_front = sys.part.select(
                lambda p: p.vs_relative[0] == self.attributes.realz_indices[x])
            slice_master_back = sys.part.select(
                lambda p: p.vs_relative[0] == self.attributes.realz_indices[x+1])
            slice_front = [x.id for x in slice_master_front if x.type ==
                           2 and x.id in self.attributes.corner_indices]
            slice_back = [x.id for x in slice_master_back if x.type ==
                          2 and x.id in self.attributes.corner_indices]
            sys.part[slice_front[1]].add_bond((self.attributes.bond, slice_back[0]))
            sys.part[slice_front[3]].add_bond((self.attributes.bond, slice_back[2]))
            sys.part[slice_front[5]].add_bond((self.attributes.bond, slice_back[4]))
            sys.part[slice_front[7]].add_bond((self.attributes.bond, slice_back[6]))

    def bond_edges(self):
        self.attributes.bond = espressomd.interactions.FeneBond(
            k=self.attributes.fene_k, d_r_max=self.attributes.fene_r_max, r_0=self.attributes.fene_r0)
        sys.bonded_inter.add(self.attributes.bond)
        for x in range(len(self.attributes.realz_indices[:-1])):
            slice_master_front = sys.part.select(
                lambda p: p.vs_relative[0] == self.attributes.realz_indices[x])
            slice_master_back = sys.part.select(
                lambda p: p.vs_relative[0] == self.attributes.realz_indices[x+1])
            slice_front = [x.id for x in slice_master_front if x.type ==
                           2 and x.id in self.attributes.edge_indices]
            slice_back = [x.id for x in slice_master_back if x.type ==
                          2 and x.id in self.attributes.edge_indices]
            sys.part[slice_front[1]].add_bond((self.attributes.bond, slice_back[0]))
            sys.part[slice_front[3]].add_bond((self.attributes.bond, slice_back[2]))
            sys.part[slice_front[5]].add_bond((self.attributes.bond, slice_back[4]))
            sys.part[slice_front[7]].add_bond((self.attributes.bond, slice_back[6]))

    def bond_relalz(self):
        self.attributes.bond = espressomd.interactions.FeneBond(
            k=self.attributes.fene_k, d_r_max=self.attributes.fene_r_max, r_0=self.attributes.fene_r0)
        sys.bonded_inter.add(self.attributes.bond)
        for x in range(len(self.attributes.realz_indices)-1):
            sys.part[self.attributes.realz_indices[x]].add_bond(
                (self.attributes.bond, self.attributes.realz_indices[x+1]))


class Simulation():

    def __init__(self, core):
        self.attributes = core

    def initialise_filament(self):
        self.filament = Filament()

    def get_sys(self):
        return self.__core

    @system_check
    def set_sys(self, core):
        self.__core = core
        globals()['sys'] = espressomd.System(box_l=core.box_l)
        seed = core.seed
        np.random.seed(seed=seed)
        sys.periodicity = core.periodicity
        sys.time_step = core.time_step
        sys.cell_system.skin = core.skin
        sys.min_global_cut = core.min_global_cut
        sys.cell_system.set_domain_decomposition
        sys.thermostat.set_langevin(
            kT=core.temperature, gamma=core.gamma, seed=seed)
        sys.virtual_sites = VirtualSitesRelative()
        # sys.seed = core.seed
        # logging.info(sys.seed)
        # sys.set_random_state_PRNG()

    attributes = property(get_sys, set_sys, None, None)

    def __getattr__(self, attrname):
        logging.info('Rerouted to Filament class, method: ' + attrname)
        return getattr(self.filament, attrname)

    @staticmethod
    def set_H_ext(H=(0, 0, 1.)):
        for x in sys.constraints:
            sys.constraints.remove(x)
            logging.info('Removed old H')

        logging.info('External field set: '+str(H))
        ExtH = espressomd.constraints.HomogeneousMagneticField(H=list(H))
        sys.constraints.add(ExtH)
        for x in sys.constraints:
            logging.info(x)

    @staticmethod
    def get_H_ext():
        for x in sys.constraints:
            logging.info(str(x.H))
            return x.H

    def init_magnetic_inter(self):
        logging.info('direct summation magnetic interactions initiated')
        dds = DipolarDirectSumCpu(prefactor=1)
        sys.actors.add(dds)

    def init_lb_GPU(self, timestep):
        logging.info('CPU LB method is beeing initiated')
        sys.thermostat.turn_off()
        sys.part[:].v = (0, 0, 0)
        lbf = espressomd.lb.LBFluidWalberla(agrid=1, dens=1, visc=args.eta, tau=timestep)
        print(sys.actors.active_actors)
        if len(sys.actors.active_actors) == 2:
            sys.actors.remove(sys.actors.active_actors[-1])
        sys.time_step = timestep
        sys.actors.add(lbf)
        # gamma_MD = 1/(1/(6*np.pi*0.25*(args.eta))-1/((args.eta)*25))
        gamma_MD = 3
        print('gamma_MD: '+str(gamma_MD))
        sys.thermostat.set_lb(LB_fluid=lbf, gamma=gamma_MD, seed=core.seed)
        logging.info(lbf.get_params())
        logging.info('GPU LB method is set with the params above.')
        return lbf

    def create_flow_channel(self, slip_vel=(0, 0, 0)):
        logging.info("Setup LB boundaries.")
        top_wall = espressomd.shapes.Wall(normal=[1, 0, 0], dist=1)
        bottom_wall = espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(140 - 1))

        top_boundary = espressomd.lbboundaries.LBBoundary(shape=top_wall, velocity=slip_vel)
        bottom_boundary = espressomd.lbboundaries.LBBoundary(shape=bottom_wall)

        sys.lbboundaries.add(top_boundary)
        sys.lbboundaries.add(bottom_boundary)

    def writevtk(self, path):
        pos = [x.pos for x in sys.part[:]]
        directors = [x.director for x in sys.part[:]]
        with open(path, 'w') as vtk:
            vtk.write("# vtk DataFile Version 2.0\n")
            vtk.write("particles\n")
            vtk.write("ASCII\n")
            vtk.write("DATASET UNSTRUCTURED_GRID\n")
            vtk.write("POINTS {} floats\n".format(len(pos)))
            for i in range(len(pos)):
                vtk.write("%f %f %f\n" % (pos[i][0], pos[i][1], pos[i][2]))

            vtk.write("POINT_DATA {}\n".format(len(pos)))
            vtk.write("SCALARS dipoles float 3\n")
            vtk.write("LOOKUP_TABLE default\n")
            for i in range(len(pos)):
                vtk.write("%f %f %f\n" % (directors[i][0], directors[i][1], directors[i][2]))
        return

    def run_anneling(self, steps, sample):
        sys.integrator.run(steps=0)
        for x in range(steps):
            logging.info('integrating step: '+str(x+1)+'/'+str(steps))
            sys.integrator.run(steps=sample)

    def filter_part_properties(self):
        # allways take the first particle with index 0
        prop_dict = dir(sys.part[1])
        stuffz = [x for x in prop_dict if (isinstance(getattr(sys.part[1], x), Iterable)
                                           or isinstance(getattr(sys.part[1], x), numbers.Number))]
        return stuffz

    def run_measurment(self, steps, sample):
        assert not self.filament.attributes.magnetize, "you are running with the non-spara run funct !!!!!"
        sys.integrator.run(steps=0)
        dict_of_god = dict()
        prop_dict = self.filter_part_properties()
        for timestep in range(steps):
            logging.info('integrating step: '+str(timestep))
            sys.integrator.run(steps=sample)
            dict_of_god['timestep_%s' % timestep] = {
                'pid_%s' % x.id: {y: getattr(x, y) for y in prop_dict} for x in sys.part[
                    self.filament.attributes.realz_indices]}
            vtf.writevcf(sys, measurments_trj)
            measurments_trj.flush()
        pickle.dump(dict_of_god, gzip.open(path_to_custom_data, 'wb'), pickle.HIGHEST_PROTOCOL)

    def run_measurment_spara(self, steps, sample):
        assert self.filament.attributes.magnetize, "you are running with the spara run funct but you shouldnt!!!!!"
        sys.integrator.run(steps=0)
        dict_of_god = dict()
        prop_dict = self.filter_part_properties()
        for timestep in range(steps):
            logging.info('integrating step: '+str(timestep))
            for bla in range(sample):
                self.filament.magnetize()
                sys.integrator.run(0, recalc_forces=True)
                sys.integrator.run(steps=1, reuse_forces=True)
            dict_of_god['timestep_%s' % timestep] = {'pid_%s' % x.id: {y: getattr(
                x, y) for y in prop_dict} for x in sys.part.select(
                    lambda p: p.type == self.filament.attributes.type_to_be_magnetized)}
            vtf.writevcf(sys, measurments_trj)
            measurments_trj.flush()
        pickle.dump(dict_of_god, gzip.open(path_to_custom_data, 'wb'), pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def init_checkpoint(lb_handle, path):
        checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path=path)
        checkpoint.register("sys.part")
        return checkpoint

    @staticmethod
    def load_checkpoint(lbb_handle, path):
        logging.info('loading checkpoint and fluid state')
        checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path=path)
        try:
            checkpoint.load()
            input = gzip.GzipFile(path_fluid_compressed, 'rb')
            s = input.read()
            input.close()
            output = open(path_fluid_sacrifice, 'wb')
            output.write(s)
            output.close()
            lbb_handle.load_checkpoint(path=path_fluid_sacrifice, binary=True)
            os.remove(path_fluid_sacrifice)
        except Exception:
            exc_type, value, traceback = sysos.exc_info()
            logging.info("Failed with exception [%s,%s ,%s]" %
                         (exc_type, value, traceback))
            sys.part.clear()
            logging.info('Retrying to load')
            checkpoint.load()
        return checkpoint

    @staticmethod
    def load_fluid(lbb_handle, path):
        logging.info('loading fluid state')
        input = gzip.GzipFile(path_fluid_compressed, 'rb')
        s = input.read()
        input.close()
        output = open(path_fluid_sacrifice, 'wb')
        output.write(s)
        output.close()
        lbb_handle.load_checkpoint(path=path_fluid_sacrifice, binary=True)
        os.remove(path_fluid_sacrifice)


# %%
dir_path = os.getcwd()
path_to_custom_data = os.path.join(
    '/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/', 'custom_dict_cube.p.gz')
measurments_trj = open(
    os.path.join('/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/', 'trajectory_measured_cube.vtf'), mode='w+t')

path_fluid_save = os.path.join(
    '/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/', 'fluid_field.txt')
path_fluid_compressed = os.path.join(
    '/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/', 'fluid_compressed.txt.gz')
path_fluid_sacrifice = os.path.join(
    '/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/', 'you_cant_see_me.txt')

core_inst = core()
sim_inst = Simulation(core_inst)
sim_inst.initialise_filament()
sim_inst.set_filament()
sim_inst.bond_cornerz()
sim_inst.bond_edges()
# sim_inst.init_magnetic_inter()
sim_inst.set_steric()
sim_inst.excommunicado()

vtf.writevsf(sys, measurments_trj)
vtf.writevcf(sys, measurments_trj)

sim_inst.set_rot_intertia_from_units()

lbb = sim_inst.init_lb_GPU(timestep=core.time_step)
sim_inst.create_flow_channel()

# # //////////////////////////////////////////////////////////////////////////////////////////////
# sys.part.clear()
# lokacija = '/Users/avantgarde_jack/PyEspresso/espresso_walberla/test_scripts/'
# path_fluid_save = os.path.join(lokacija, 'fluid_field.txt')
# path_fluid_compressed = os.path.join(lokacija, 'fluid_compressed.txt.gz')
# path_fluid_sacrifice = os.path.join(lokacija, 'you_cant_see_me.txt')
# # sim_inst.load_fluid(lbb, path=lokacija)
# checkpoint = sim_inst.load_checkpoint(lbb, path=lokacija)
# sys.integrator.run(steps=0, reuse_forces=True)
# logging.info("Loaded LB fluid state!")
# # /////////////////////////////////////////////////////////////////////////////////////////////

print(sim_inst.filament.attributes.init_pos)
today = datetime.now()
print("Today's date:", today)

logging.info("Warming up the system with LB fluid.")
for x in range(100):
    sys.integrator.run(600)
    logging.info("converging slowly.")
    vtf.writevcf(sys, measurments_trj)
logging.info("Warming up the system with LB fluid finished.")

# if args.h_field:
#     sim_inst.set_H_ext(H=(0, 0, args.h_field))
sys.lbboundaries.clear()
sim_inst.create_flow_channel(slip_vel=args.slip_v)

sim_inst.run_measurment(500, 5000)

# logging.info("Running fluid velocity measurments")
# for x in range(10):
#     sys.integrator.run(60000)
#     logging.info("getting there.")
#     lbb.print_vtk_velocity(os.path.join(os.path.join(
#         '/home/lv70806/', dir_path[14:]), 'vtk_velo_' + str(x) + ".vtk"))
# logging.info("Warming up the system with LB fluid finished.")
# sim_inst.integrator_ferro()

# # //////////////////////////////////////////////////////////////////////////////////////////////
# checkpoint=sim_inst.init_checkpoint(lbb,path=lokacija)
# checkpoint.save()
# lbb.save_checkpoint(path=path_fluid_save, binary=True)
# with open(path_fluid_save, 'rb') as f_in, gzip.open(path_fluid_compressed, 'wb') as f_out:
#     shutil.copyfileobj(f_in, f_out)
# os.remove(path_fluid_save)
# # /////////////////////////////////////////////////////////////////////////////////////////////

today = datetime.now()
print("Today's date:", today)
logging.info('done with simulation instance!')
