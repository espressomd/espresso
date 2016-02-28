#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

import h5py
import sys
import numpy as np
import espressomd._system as es
from espressomd import analyze
include "myconfig.pxi"


class h5md(object):

    def __init__(self, filename, system):
        self.system = system
        self.filename = filename
        self.h5_file = h5py.File(filename, "a")
        self.write_to_h5 = self.write_to_h5(self)
        self.read_from_h5 = self.read_from_h5(self)

    def close(self):
        self.h5_file.close()

    def WriteValueUserdefined(
            self, value, h5_datasetgroup, h5_index, h5_maxshape, h5_Dtype, chunk=True):
        # value:             Value to be written
        # h5_datasetgroup:    Path and name of h5 dataset
        # h5_index:          Index of the to be written value
        # h5_maxshape:       Maximum shape size
        # h5_Dtype:          Datatype of h5 dataset
        # chunk:             Chunk size of h5 dataset

        # DATASET SHAPE
        h5_shape = ()
        for i in range(len(h5_index)):
            h5_shape = h5_shape + (h5_index[i] + 1,)

        # CREATE OR OPEN DATASET
        try:
            self.dataset = self.h5_file.create_dataset(
                h5_datasetgroup, h5_shape, maxshape=h5_maxshape, dtype=h5_Dtype, chunks=chunk)
        except:
            self.dataset = self.h5_file[h5_datasetgroup]

        # RESIZE
        for i in range(len(h5_index)):
            if(self.dataset.shape[i] < h5_shape[i]):
                self.dataset.resize(h5_shape)
                break

        # ASSIGN VALUE
        self.dataset[h5_index] = value

    def WriteValueEspresso(self, timestep, particle_id, value, h5_datasetgroup,
                           h5_shape, h5_maxshape, h5_Dtype, chunk=True, case="", feature=1):
        # timestep:          ESPResSo time step
        # particle_id:       ESPResSo particle id
        # value:             Value to be written
        # h5_datasetgroup:    Path and name of h5 dataset
        # h5_shape:          Shape of h5 dataset
        # h5_maxshape:       Maximum shape size
        # h5_Dtype:          Datatype of h5 dataset
        # chunk:             Chunk size of h5 dataset
        # case:              Distinguish different write cases to avoid multiple WriteValue-functions
        # feature:           Required ESPResSo features

        # CREATE OR OPEN DATASET
        if feature != 1:
            print "ERROR H5: Some necessary ESPResSo features for values used in h5-file are not activated"
            sys.exit()
        try:
            self.dataset = self.h5_file.create_dataset(
                h5_datasetgroup, h5_shape, maxshape=h5_maxshape, dtype=h5_Dtype, chunks=chunk)
        except:
            self.dataset = self.h5_file[h5_datasetgroup]

        # WRITE CASES:
        #time and step
        if case == 'time_dependent':
            # Resize dataset if new time step added and assign value
            maxshape_timestep = max((self.dataset.shape[0], timestep + 1))
            if(self.dataset.shape[0] < timestep + 1):
                self.dataset.resize((maxshape_timestep, 1))
            self.dataset[timestep] = value

        # particle properties
        if case == 'particle_time_independent':
            maxshape_n_part = max((self.dataset.shape[0], self.system.n_part))
            if(self.dataset.shape[0] < self.system.n_part):
                self.dataset.resize((maxshape_n_part, self.dataset.shape[1]))
            self.dataset[particle_id] = value
        if case == 'particle_time_dependent':
            maxshape_timestep = max((self.dataset.shape[0], timestep + 1))
            maxshape_n_part = max((self.dataset.shape[1], self.system.n_part))
            if(self.dataset.shape[0] < timestep + 1 or self.dataset.shape[1] < self.system.n_part):
                self.dataset.resize(
                    (maxshape_timestep, maxshape_n_part, self.dataset.shape[2]))
            self.dataset[timestep, particle_id] = value
                                                
        # particle properties bond
        if case == 'particle_bond_time_independent':
            maxshape_particle_id = max(
                (self.dataset.shape[0], particle_id + 1))
            if(self.dataset.shape[0] < particle_id + 1):
                self.dataset.resize((maxshape_particle_id,))
            self.dataset[particle_id] = value
        if case == 'particle_bond_time_dependent':
            maxshape_timestep = max((self.dataset.shape[0], timestep + 1))
            maxshape_particle_id = max(
                (self.dataset.shape[1], particle_id + 1))
            if(self.dataset.shape[0] < timestep + 1 or self.dataset.shape[1] < particle_id + 1):
                self.dataset.resize(
                    (maxshape_timestep, maxshape_particle_id, 1))
            self.dataset[timestep, particle_id] = value

        # observables
        if case == 'observable_time_independent':
            value_temp = []
            try:  # value is array
                value_length = len(value)
                value_temp = value
            except:  # value is scalar
                value_length = 1
                value_temp.append(value)
            self.dataset.resize((value_length,))
            for i in range(0, value_length):
                self.dataset[i] = value_temp[i]
        if case == 'observable_time_dependent':
            value_temp = []
            try:  # value is array
                try:  # Transform matrix to array
                    value = value.ravel()
                except:
                    pass
                value_length = len(value)
                value_temp = value
            except:  # value is scalar
                value_length = 1
                value_temp.append(value)
            maxshape_timestep = max((self.dataset.shape[0], timestep + 1))
            if(self.dataset.shape[0] < timestep + 1):
                self.dataset.resize((maxshape_timestep, value_length))
            self.dataset[timestep] = value_temp

        # box
        if case == 'box_edges_time_independent':
            self.dataset[0, 0] = self.system.box_l[0]
            self.dataset[1, 1] = self.system.box_l[1]
            self.dataset[2, 2] = self.system.box_l[2]
        if case == 'box_edges_time_dependent':
            if(self.dataset.shape[0] <= timestep + 1):
                self.dataset.resize((timestep + 1, 3, 3))
            self.dataset[timestep, 0, 0] = self.system.box_l[0]
            self.dataset[timestep, 1, 1] = self.system.box_l[1]
            self.dataset[timestep, 2, 2] = self.system.box_l[2]
        if case == 'box_dimension_time_independent':
            self.dataset[0] = value
        if case == 'box_dimension_time_dependent':
            if(self.dataset.shape[0] <= timestep + 1):
                self.dataset.resize((timestep + 1, 1))
            self.dataset[timestep, 0] = value
        if case == 'box_boundary_time_independent':
            if(value[0] == 0):
                self.dataset[0] = "none"
            else:
                self.dataset[0] = "periodic"
            if(value[1] == 0):
                self.dataset[1] = "none"
            else:
                self.dataset[1] = "periodic"
            if(value[2] == 0):
                self.dataset[2] = "none"
            else:
                self.dataset[2] = "periodic"
        if case == 'box_boundary_time_dependent':
            if(self.dataset.shape[0] <= timestep + 1):
                self.dataset.resize((timestep + 1, 3))
            if(value[0] == 0):
                self.dataset[timestep, 0] = "none"
            else:
                self.dataset[timestep, 0] = "periodic"
            if(value[1] == 0):
                self.dataset[timestep, 1] = "none"
            else:
                self.dataset[timestep, 1] = "periodic"
            if(value[2] == 0):
                self.dataset[timestep, 2] = "none"
            else:
                self.dataset[timestep, 2] = "periodic"


#----------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------WRITE CLASS------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
# WRITE CLASS
    class write_to_h5(object):

        def __init__(self, h5md):
            self.h5md = h5md

    # USERDEFINED
        # user defined dataset
        def userdefined(
                self, value, groupname, datasetname, datatype, dataset_index, chunk=True):
            self.h5md.WriteValueUserdefined(
                value, groupname + datasetname, dataset_index, (None), datatype, chunk)

    # TIME
        # time
        def time(self, timestep=-1, groupname="particles/atoms/Time/",
                 datasetname="time", chunk=True):
            self.h5md.WriteValueEspresso(
                timestep, -1, self.h5md.system.time, groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'time_dependent', 1)
        # time step

        def time_step(
                self, timestep=-1, groupname="particles/atoms/Step/", datasetname="step", chunk=True):
            self.h5md.WriteValueEspresso(
                timestep, -1, timestep, groupname + datasetname, (1, 1), (None, 1), 'int32', chunk, 'time_dependent', 1)

    # PARTICLE PROPERTIES
        # position
        def pos(self, timestep=-1, groupname="particles/atoms/position/",
                datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].pos, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].pos, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', 1)
                    self.time(timestep, groupname, "time", chunk)
                    self.time_step(timestep, groupname, "step", chunk)
                    
        # velocity
        def v(self, timestep=-1, groupname="particles/atoms/velocity/",
              datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].v, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].v, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', 1)
        
        # force
        def f(self, timestep=-1, groupname="particles/atoms/force/",
              datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].f, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].f, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', 1)
        
        # bonds
        def bonds(self, timestep=-1, groupname="particles/atoms/",
                  datasetname="value", chunk=True):
            self.bond_from(
                timestep, groupname + "bond_from/", datasetname, chunk)
            self.bond_to(timestep, groupname + "bond_to/", datasetname, chunk)

        def bond_from(
                self, timestep=-1, groupname="particles/atoms/bond_from/", datasetname="value", chunk=True):
            index = 0
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    for bond in self.h5md.system.part[i].bonds:
                        for partner in range(1, len(bond)):
                            self.h5md.WriteValueEspresso(-1, index, i, groupname + datasetname, (
                                1,), (None,), 'int32', chunk, 'particle_bond_time_independent', 1)
                            index += 1
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    for bond in self.h5md.system.part[i].bonds:
                        for partner in range(1, len(bond)):
                            self.h5md.WriteValueEspresso(
                                timestep, index, i, groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_bond_time_dependent', 1)
                            index += 1

        def bond_to(self, timestep=-1, groupname="particles/atoms/bond_to/",
                    datasetname="value", chunk=True):
            index = 0
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    for bond in self.h5md.system.part[i].bonds:
                        for partner in range(1, len(bond)):
                            self.h5md.WriteValueEspresso(-1, index, bond[partner], groupname + datasetname, (
                                1,), (None,), 'int32', chunk, 'particle_bond_time_independent', 1)
                            index += 1
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    for bond in self.h5md.system.part[i].bonds:
                        for partner in range(1, len(bond)):
                            self.h5md.WriteValueEspresso(timestep, index, bond[
                                                         partner], groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_bond_time_dependent', 1)
                            index += 1
        
        # type / species
        def type(self, timestep=-1, groupname="particles/atoms/type/",
                 datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].type, groupname + datasetname, (1, 1), (None, 1), 'int32', chunk, 'particle_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].type, groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_time_dependent', 1)
                            
        # id
        def id(self, timestep=-1, groupname="particles/atoms/id/",
               datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].id, groupname + datasetname, (1, 1), (None, 1), 'int32', chunk, 'particle_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].id, groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_time_dependent', 1)
        
        # mass
        def mass(self, timestep=-1, groupname="particles/atoms/mass/",
                 datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].mass, groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'particle_time_independent', MASS)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].mass, groupname + datasetname, (1, 1, 1), (None, None, 1), 'f8', chunk, 'particle_time_dependent', MASS)
        
        # omega_lab
        def omega_lab(
                self, timestep=-1, groupname="particles/atoms/omega_lab/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].omega_lab, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', ROTATION)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].omega_lab, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', ROTATION)
        
        # rinertia
        def rinertia(self, timestep=-1, groupname="particles/atoms/rinertia/",
                     datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].rinertia, groupname + datasetname, (
                        1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', ROTATIONAL_INERTIA)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].rinertia, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', ROTATIONAL_INERTIA)
        
        # omega_body
        def omega_body(
                self, timestep=-1, groupname="particles/atoms/omega_body/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].omega_body, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', ROTATION)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].omega_body, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', ROTATION)
        
        # torque_lab
        def torque_lab(
                self, timestep=-1, groupname="particles/atoms/torque_lab/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].torque_lab, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', ROTATION)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].torque_lab, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', ROTATION)
        
        # quat
        def quat(self, timestep=-1, groupname="particles/atoms/quat/",
                 datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].quat, groupname + datasetname, (1, 4), (None, 4), 'f8', chunk, 'particle_time_independent', ROTATION)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].quat, groupname + datasetname, (1, 1, 4), (None, None, 4), 'f8', chunk, 'particle_time_dependent', ROTATION)
        
        # charge
        def q(self, timestep=-1, groupname="particles/atoms/charge/",
              datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].q, groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'particle_time_independent', ELECTROSTATICS)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].q, groupname + datasetname, (1, 1, 1), (None, None, 1), 'f8', chunk, 'particle_time_dependent', ELECTROSTATICS)
        
        # virtual
        def virtual(self, timestep=-1, groupname="particles/atoms/virtual/",
                    datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].virtual, groupname + datasetname, (
                        1, 1), (None, 1), 'int32', chunk, 'particle_time_independent', VIRTUAL_SITES_COM)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].virtual, groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_time_dependent', VIRTUAL_SITES_COM)
        
        # vs_relative
        def vs_relative(
                self, timestep=-1, groupname="particles/atoms/vs_relative/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].vs_relative, groupname + datasetname, (
                        1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', VIRTUAL_SITES_RELATIVE)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].vs_relative, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', VIRTUAL_SITES_RELATIVE)
        
        # dipole
        def dip(self, timestep=-1, groupname="particles/atoms/dipole/",
                datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].dip, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', DIPOLES)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].dip, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', DIPOLES)
        
        # dipole_magnitude
        def dipm(self, timestep=-1, groupname="particles/atoms/dipole_magnitude/",
                 datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].dipm, groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'particle_time_independent', DIPOLES)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].dipm, groupname + datasetname, (1, 1, 1), (None, None, 1), 'f8', chunk, 'particle_time_dependent', DIPOLES)
        
        # external force
        def ext_force(
                self, timestep=-1, groupname="particles/atoms/ext_force/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].ext_force, groupname + datasetname, (
                        1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', EXTERNAL_FORCES)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].ext_force, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', EXTERNAL_FORCES)
        
        # external force particle fix
        def fix(self, timestep=-1, groupname="particles/atoms/ext_force_fix/",
                datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].fix, groupname + datasetname, (
                        1, 3), (None, 3), 'int32', chunk, 'particle_time_independent', EXTERNAL_FORCES)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].fix, groupname + datasetname, (1, 1, 3), (None, None, 3), 'int32', chunk, 'particle_time_dependent', EXTERNAL_FORCES)
        
        # external torque
        def ext_torque(
                self, timestep=-1, groupname="particles/atoms/ext_torque/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[
                                                 i].ext_torque, groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'particle_time_independent', ROTATION)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].ext_torque, groupname + datasetname, (1, 1, 3), (None, None, 3), 'f8', chunk, 'particle_time_dependent', ROTATION)
        
        # gamma
        def gamma(self, timestep=-1, groupname="particles/atoms/gamma/",
                  datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].gamma, groupname + datasetname, (
                        1, 1), (None, 1), 'f8', chunk, 'particle_time_independent', LANGEVIN_PER_PARTICLE)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].gamma, groupname + datasetname, (1, 1, 1), (None, None, 1), 'f8', chunk, 'particle_time_dependent', LANGEVIN_PER_PARTICLE)
        
        # temperature
        def temp(self, timestep=-1, groupname="particles/atoms/temp/",
                 datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].temp, groupname + datasetname, (
                        1, 1), (None, 1), 'f8', chunk, 'particle_time_independent', LANGEVIN_PER_PARTICLE)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].temp, groupname + datasetname, (1, 1, 1), (None, None, 1), 'f8', chunk, 'particle_time_dependent', LANGEVIN_PER_PARTICLE)
        
        # rotation
        def rotation(self, timestep=-1, groupname="particles/atoms/rotation/",
                     datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(-1, i, self.h5md.system.part[i].rotation, groupname + datasetname, (
                        1, 1), (None, 1), 'int32', chunk, 'particle_time_independent', ROTATION_PER_PARTICLE)
            if timestep >= 0:  # Time dependent
                for i in range(0, self.h5md.system.n_part):
                    self.h5md.WriteValueEspresso(timestep, i, self.h5md.system.part[
                                                 i].rotation, groupname + datasetname, (1, 1, 1), (None, None, 1), 'int32', chunk, 'particle_time_dependent', ROTATION_PER_PARTICLE)

    # OBSERVABLES
        # mindist
        def mindist(self, timestep=-1, p1='default', p2='default',
                    groupname="observables/mindist/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.mindist(self.h5md.system, p1, p2),
                                             groupname + datasetname, (1,), (1,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.mindist(self.h5md.system, p1, p2),
                                             groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'observable_time_dependent', 1)
        
        # distto
        def distto(self, timestep=-1, id=None, pos=None,
                   groupname="observables/distto/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.distto(self.h5md.system, id, pos),
                                             groupname + datasetname, (1,), (1,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.distto(self.h5md.system, id, pos),
                                             groupname + datasetname, (1, 1), (None, 1), 'f8', chunk, 'observable_time_dependent', 1)
        
        # analyze_linear_momentum
        def analyze_linear_momentum(self, timestep=-1, include_particles=True, include_lbfluid=True,
                                    groupname="observables/linear_momentum/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.analyze_linear_momentum(self.h5md.system, include_particles,
                                                                                           include_lbfluid), groupname + datasetname, (3,), (3,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.analyze_linear_momentum(self.h5md.system, include_particles,
                                                                                           include_lbfluid), groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'observable_time_dependent', 1)
        
        # nbhood
        def nbhood(self, timestep=-1, pos=None, r_catch=None, plane='3d',
                   groupname="observables/nbhood/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.nbhood(self.h5md.system, pos, r_catch, plane),
                                             groupname + datasetname, (4,), (None,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.nbhood(self.h5md.system, pos, r_catch, plane),
                                             groupname + datasetname, (1, 4), (None, None), 'f8', chunk, 'observable_time_dependent', 1)
        
        # cylindrical_average
        def cylindrical_average(self, timestep=-1, center=None, direction=None, length=None, radius=None, bins_axial=None,
                                bins_radial=None, types=[-1], groupname="observables/cylindrical_average/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(analyze.cylindrical_average(
                        self.h5md.system, center, direction, length, radius, bins_axial, bins_radial, types).shape[0]):
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.cylindrical_average(self.h5md.system, center, direction, length, radius, bins_axial, bins_radial, types)[
                                                 i], groupname + datasetname + "_" + str(i), (8,), (None,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(analyze.cylindrical_average(
                        self.h5md.system, center, direction, length, radius, bins_axial, bins_radial, types).shape[0]):
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.cylindrical_average(self.h5md.system, center, direction, length, radius, bins_axial, bins_radial, types)[
                                                 i], groupname + datasetname + "_" + str(i), (1, 8), (None, None), 'f8', chunk, 'observable_time_dependent', 1)
        
        # pressure
        def pressure(self, timestep=-1, v_comp=False,
                     groupname="observables/pressure/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for key in analyze.pressure(self.h5md.system, v_comp).keys():
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.pressure(self.h5md.system, v_comp)[
                                                 key], groupname + datasetname + "_" + str(key), (1,), (None,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for key in analyze.pressure(self.h5md.system, v_comp).keys():
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.pressure(self.h5md.system, v_comp)[
                                                 key], groupname + datasetname + "_" + str(key), (1, 1), (None, 1), 'f8', chunk, 'observable_time_dependent', 1)
        
        # stress_tensor
        def stress_tensor(self, timestep=-1, v_comp=False,
                          groupname="observables/stress_tensor/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for key in analyze.stress_tensor(
                        self.h5md.system, v_comp).keys():
                    try:  # value is array
                        if(len(analyze.stress_tensor(self.h5md.system, v_comp)[key].shape) > 1):
                            for i in range(
                                    analyze.stress_tensor(self.h5md.system, v_comp)[key].shape[0]):
                                self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[key][
                                                             i], groupname + datasetname + "_" + str(key) + "_" + str(i), (9,), (9,), 'f8', chunk, 'observable_time_independent', 1)
                        else:
                            self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[
                                                         key], groupname + datasetname + "_" + str(key), (9,), (9,), 'f8', chunk, 'observable_time_independent', 1)
                    except:  # value is scalar
                        self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[
                                                     key], groupname + datasetname + "_" + str(key), (9,), (9,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for key in analyze.stress_tensor(
                        self.h5md.system, v_comp).keys():
                    try:  # value is array
                        if(len(analyze.stress_tensor(self.h5md.system, v_comp)[key].shape) > 1):
                            for i in range(
                                    analyze.stress_tensor(self.h5md.system, v_comp)[key].shape[0]):
                                self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[key][
                                                             i], groupname + datasetname + "_" + str(key) + "_" + str(i), (1, 9), (None, 9), 'f8', chunk, 'observable_time_dependent', 1)
                        else:
                            self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[
                                                         key], groupname + datasetname + "_" + str(key), (1, 9), (None, 9), 'f8', chunk, 'observable_time_dependent', 1)
                    except:  # value is scalar
                        self.h5md.WriteValueEspresso(timestep, -1, analyze.stress_tensor(self.h5md.system, v_comp)[
                                                     key], groupname + datasetname + "_" + str(key), (1, 9), (None, 9), 'f8', chunk, 'observable_time_dependent', 1)
        
        # local_stress_tensor
        def local_stress_tensor(self, timestep=-1, system=None, periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(
                1.0, 1.0, 1.0), bins=(1, 1, 1), groupname="observables/pressure/", datasetname="value", chunk=True):
            print("ERROR H5: WRITE local_stress_tensor not implemented yet")
        
        # energy
        def analyze_energy(self, timestep=-1, etype='all', id1='default', id2='default',
                           groupname="observables/energy/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for key in analyze.energy(self.h5md.system).keys():
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.energy(self.h5md.system, etype, id1, id2)[
                                                 key], groupname + datasetname + "_" + str(key), (1,), (None,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for key in analyze.energy(self.h5md.system).keys():
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.energy(self.h5md.system, etype, id1, id2)[
                                                 key], groupname + datasetname + "_" + str(key), (1, 1), (None, None), 'f8', chunk, 'observable_time_dependent', 1)
        
        # calc_re
        def calc_re(self, timestep=-1, chain_start=None, number_of_chains=None,
                    chain_length=None, groupname="observables/calc_re/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_re(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (3,), (3,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_re(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'observable_time_dependent', 1)
        
        # calc_rg
        def calc_rg(self, timestep=-1, chain_start=None, number_of_chains=None,
                    chain_length=None, groupname="observables/calc_rg/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_rg(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (3,), (3,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_rg(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'observable_time_dependent', 1)
        
        # calc_rh
        def calc_rh(self, timestep=-1, chain_start=None, number_of_chains=None,
                    chain_length=None, groupname="observables/calc_rh/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_rh(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (3,), (3,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, analyze.calc_rh(self.h5md.system, chain_start, number_of_chains,
                                                                           chain_length), groupname + datasetname, (1, 3), (None, 3), 'f8', chunk, 'observable_time_dependent', 1)
        
        # check_topology
        def check_topology(self, timestep=-1, system=None, chain_start=None, number_of_chains=None,
                           chain_length=None, groupname="observables/check_topology/", datasetname="value"):
            print("ERROR H5: check_topology not implemented yet")
        
        # structure_factor
        def structure_factor(self, timestep=-1, sf_type='default', sf_order='default',
                             groupname="observables/structure_factor/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                for i in range(
                        len(analyze.structure_factor(self.h5md.system, sf_type, sf_order))):
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.structure_factor(self.h5md.system, sf_type, sf_order)[
                                                 i], groupname + datasetname + "_" + str(i), (2,), (2,), 'f8', chunk, 'observable_time_independent', 1)
            if timestep >= 0:  # Time dependent
                for i in range(
                        len(analyze.structure_factor(self.h5md.system, sf_type, sf_order))):
                    self.h5md.WriteValueEspresso(timestep, -1, analyze.structure_factor(self.h5md.system, sf_type, sf_order)[
                                                 i], groupname + datasetname + "_" + str(i), (1, 2), (None, 2), 'f8', chunk, 'observable_time_dependent', 1)
    
    # BOX
        def box_edges(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="value", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(-1, -1, self.h5md.system.box_l, groupname +
                                             datasetname, (3, 3), (3, 3), 'f8', chunk, 'box_edges_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, self.h5md.system.box_l, groupname +
                                             datasetname, (1, 3, 3), (None, 3, 3), 'f8', chunk, 'box_edges_time_dependent', 1)

        def box_boundary(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="boundary", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(
                    timestep, -1, self.h5md.system.periodicity, groupname + datasetname, (3,), (3,), 'S30', chunk, 'box_boundary_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(timestep, -1, self.h5md.system.periodicity, groupname +
                                             datasetname, (1, 3), (None, 3), 'S30', chunk, 'box_boundary_time_dependent', 1)

        def box_dimension(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="dimension", chunk=True):
            if timestep == -1:  # Time independent
                self.h5md.WriteValueEspresso(
                    timestep, -1, 3, groupname + datasetname, (1,), (1,), 'int32', chunk, 'box_dimension_time_independent', 1)
            if timestep >= 0:  # Time dependent
                self.h5md.WriteValueEspresso(
                    timestep, -1, 3, groupname + datasetname, (1, 1), (None, 1), 'int32', chunk, 'box_dimension_time_dependent', 1)
    
    # VMD
        def VMD(self, datasetname, value=-1,groupname="parameters/vmd_structure/", chunk=True):
            if(datasetname == 'species'):
                Value=[]
                for i in range(0, self.h5md.system.n_part):
                    Value.append(self.h5md.system.part[i].type)
                try:
                    self.dataset = self.h5md.h5_file.create_dataset("particles/atoms/species", (len(Value),), maxshape=(None), dtype='int32', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file["particles/atoms/species"]
                self.dataset.resize((len(Value),))
                for i in range(len(Value)):
                    self.dataset[i] = Value[i]
                                                            
            if(datasetname == 'charge'):
                Value=[]
                for i in range(0, self.h5md.system.n_part):
                    Value.append(self.h5md.system.part[i].q)
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(groupname + "charge", (len(Value),), maxshape=(None), dtype='f8', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "charge"]
                self.dataset.resize((len(Value),))
                for i in range(len(Value)):
                    self.dataset[i] = Value[i]
                
            if(datasetname == 'mass'):
                Value=[]
                for i in range(0, self.h5md.system.n_part):
                    Value.append(self.h5md.system.part[i].mass)
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(groupname + "mass", (len(Value),), maxshape=(None), dtype='f8', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "mass"]
                self.dataset.resize((len(Value),))
                for i in range(len(Value)):
                    self.dataset[i] = Value[i]
                
            if(datasetname == 'bond_from'):
                self.bond_from(-1, groupname, "bond_from", chunk)
                
            if(datasetname == 'bond_to'):
                self.bond_to(-1, groupname, "bond_to", chunk)

            if(datasetname == 'indexOfSpecies'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "indexOfSpecies", (len(value),), maxshape=(None), dtype='int32', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[
                        groupname + "indexOfSpecies"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'radius'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "radius", (len(value),), maxshape=(None), dtype='f8', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "radius"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'name'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "name", (len(value),), maxshape=(None), dtype='S30', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "name"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'type'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "type", (len(value),), maxshape=(None), dtype='S30', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "type"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'resname'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "resname", (len(value),), maxshape=(None), dtype='S30', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "resname"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'resid'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "resid", (len(value),), maxshape=(None), dtype='int32', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "resid"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'segid'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "segid", (len(value),), maxshape=(None), dtype='S30', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "segid"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'chain'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "chain", (len(value),), maxshape=(None), dtype='S30', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[groupname + "chain"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]
            if(datasetname == 'atomicnumber'):
                try:
                    self.dataset = self.h5md.h5_file.create_dataset(
                        groupname + "atomicnumber", (len(value),), maxshape=(None), dtype='int32', chunks=chunk)
                except:
                    self.dataset = self.h5md.h5_file[
                        groupname + "atomicnumber"]
                self.dataset.resize((len(value),))
                for i in range(len(value)):
                    self.dataset[i] = value[i]



#----------------------------------------------------------------------------------------------------------------#
#------------------------------------------------------READ CLASS------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------#
# READ CLASS
    class read_from_h5(object):

        def __init__(self, h5md):
            self.h5md = h5md

    # USERDEFINED
        def userdefined(self, groupname, datasetname, h5_index):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()

            # Read value
            return self.h5md.value_dataset[h5_index]

    # PARTICLES
        # time
        def time(self, timestep=-1, groupname="particles/atoms/Time/",
                 datasetname="time"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                self.h5md.system.time = self.h5md.value_dataset[timestep]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # time
        def time_step(
                self, timestep=-1, groupname="particles/atoms/Timestep/", datasetname="step"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                return self.h5md.value_dataset[timestep]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # Position
        def pos(self, timestep=-1, groupname="particles/atoms/position/",
                datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].pos = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].pos = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # velocity
        def v(self, timestep=-1, groupname="particles/atoms/velocity/",
              datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[i].v = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].v = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # force
        def f(self, timestep=-1, groupname="particles/atoms/force/",
              datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[i].f = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].f = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # bonds
        def bonds(
                self, timestep=-1, groupname="particles/atoms/bond_from/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset_from = self.h5md.h5_file[
                    "particles/atoms/bond_from/"][datasetname]
                self.h5md.value_dataset_to = self.h5md.h5_file[
                    "particles/atoms/bond_to/"][datasetname]
            except:
                print "ERROR H5: No " + "particles/atoms/bonds_from_to/" + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                print("ERROR H5: Reading bonds not implemented yet")
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # species
        def type(self, timestep=-1, groupname="particles/atoms/type/",
                 datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[i].type = int(self.h5md.value_dataset[i][0])
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[i].type = int(self.h5md.value_dataset[timestep, i][0])
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # id
        def id(self, timestep=-1, groupname="particles/atoms/id/",
               datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].id = int(self.h5md.value_dataset[i][0])
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].id = int(self.h5md.value_dataset[timestep, i][0])
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()

        def mass(self, timestep=-1, groupname="particles/atoms/mass/",
                 datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].mass = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].mass = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # omega_lab
        def omega_lab(
                self, timestep=-1, groupname="particles/atoms/omega_lab/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].omega_lab = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].omega_lab = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # rinertia
        def rinertia(
                self, timestep=-1, groupname="particles/atoms/rinertia/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].rinertia = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].rinertia = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # omega_body
        def omega_body(
                self, timestep=-1, groupname="particles/atoms/omega_body/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].omega_body = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].omega_body = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # torque_lab
        def torque_lab(
                self, timestep=-1, groupname="particles/atoms/torque_lab/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].torque_lab = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].torque_lab = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # quat
        def quat(self, timestep=-1, groupname="particles/atoms/quat/",
                 datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].quat = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].quat = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # charge
        def q(self, timestep=-1, groupname="particles/atoms/charge/",
              datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].q = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].q = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # virtual
        def virtual(
                self, timestep=-1, groupname="particles/atoms/virtual/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].virtual = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].virtual = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # vs_relative
        def vs_relative(
                self, timestep=-1, groupname="particles/atoms/vs_relative/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].vs_relative = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].vs_relative = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # dipole
        def dip(self, timestep=-1, groupname="particles/atoms/dipole/",
                datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].dip = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].dip = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # dipole_magnitude
        def dipm(self, timestep=-1,
                 groupname="particles/atoms/dipole_magnitude/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].dipm = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].dipm = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # external force
        def ext_force(
                self, timestep=-1, groupname="particles/atoms/ext_force/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].ext_force = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].ext_force = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # external force particle fix
        def fix(self, timestep=-1,
                groupname="particles/atoms/ext_force_fix/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].fix = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].fix = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # external torque
        def ext_torque(
                self, timestep=-1, groupname="particles/atoms/ext_torque/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].ext_torque = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].ext_torque = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # gamma
        def gamma(
                self, timestep=-1, groupname="particles/atoms/gamma/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].gamma = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].gamma = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # temperature
        def temp(self, timestep=-1, groupname="particles/atoms/temp/",
                 datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].temp = self.h5md.value_dataset[i][0]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].temp = self.h5md.value_dataset[timestep, i][0]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()
        
        # rotation
        def rotation(
                self, timestep=-1, groupname="particles/atoms/rotation/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    for i in range(self.h5md.value_dataset.shape[0]):
                        self.h5md.system.part[
                            i].rotation = self.h5md.value_dataset[i]
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    for i in range(self.h5md.value_dataset.shape[1]):
                        self.h5md.system.part[
                            i].rotation = self.h5md.value_dataset[timestep, i]
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()

    # BOX
        def box_edges(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="value"):
            # Try to open dataset
            try:
                self.h5md.value_dataset = self.h5md.h5_file[
                    groupname][datasetname]
            except:
                print "ERROR H5: No " + groupname + datasetname + " dataset in h5-file exisiting"
                sys.exit()
            # Try to read value from h5-file and write to ESPResSo
            try:
                if len(self.h5md.value_dataset.shape) == 2:  # time independent
                    box_l_temp = [-1, -1, -1]
                    for i in range(self.h5md.value_dataset.shape[0]):
                        box_l_temp[i] = self.h5md.value_dataset[i, i]
                    self.h5md.system.box_l = box_l_temp
                if len(self.h5md.value_dataset.shape) == 3:  # time dependent
                    box_l_temp = [-1, -1, -1]
                    for i in range(self.h5md.value_dataset.shape[1]):
                        box_l_temp[i] = self.h5md.value_dataset[timestep, i, i]
                    self.h5md.system.box_l = box_l_temp
            except:
                print "ERROR H5: Access to dataset " + groupname + datasetname + " or writing to ESPResSo not possible"
                sys.exit()

        def box_boundary(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="value"):
            print("ERROR H5: READ box_boundary not implemented yet")

        def box_dimension(
                self, timestep=-1, groupname="particles/atoms/box/", datasetname="value"):
            print("ERROR H5: READ box_dimension not implemented yet")
