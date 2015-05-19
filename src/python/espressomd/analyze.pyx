#
# Copyright (C) 2013,2014 The ESPResSo project
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
# For C-extern Analysis

cimport c_analyze
cimport utils
cimport particle_data
import utils
import code_info
import particle_data

#
# Minimal distance between particles
#
def mindist(system, p1 = 'default', p2 = 'default'):

  cdef IntList* set1
  cdef IntList* set2

  if p1 == 'default' and p2 == 'default':
    result = c_analyze.mindist(NULL,NULL)
  elif p1 == 'default' and not p2 == 'default':
    print 'usage: mindist([typelist],[typelist])'
    return 0
  elif not p1 == 'default' and p2 == 'default':
    print 'usage: mindist([typelist],[typelist])'
    return 0
  else:
    for i in range(len(p1)):
      if not isinstance(p1[i],int):
        print 'usage: mindist([typelist],[typelist])'
        return 0

    for i in range(len(p2)):
      if not isinstance(p2[i],int):
        print 'usage: mindist([typelist],[typelist])'
        return 0
  
    set1 = create_IntList_from_python_object(p1)
    set2 = create_IntList_from_python_object(p2)

    result = c_analyze.mindist(set1, set2)

    realloc_intlist(set1, 0)
    realloc_intlist(set2, 0)

# The following lines are probably not necessary.
#  free (set1)
#  free (set2)

  return result

#
# Distance to particle or point
#
def distto(system, id_or_pos):
  cdef double cpos[3]
  if system.n_part == 0:
    print  'no particles'
    return 'no particles'

  # check if id_or_pos is position or particle id
  if isinstance(id_or_pos,int):
    _id = id_or_pos
    _pos = particle_data.ParticleHandle(id_or_pos).pos
    for i in range(3):
      cpos[i] = _pos[i]
  else:
    for i in range(3):
      cpos[i] = id_or_pos[i]
      _id = -1
  return c_analyze.distto(cpos,_id)

# get all particles in neighborhood r_catch of pos and return their ids
# in il. plane can be used to specify the distance in the xy, xz or yz
# plane
def nbhood(self, pos=None, r_catch=None, plane = '3d'):
  cdef int planedims[3]
  cdef IntList* il=NULL
  cdef double c_pos[3]

  checkTypeOrExcept(pos, 3, float, "_pos=(float,float,float) must be passed to nbhood")
  checkTypeOrExcept(r_catch, 1, float, "r_catch=float needs to be passed to nbhood")

  #default 3d takes into account dist in x, y and z
  planedims[0] = 1
  planedims[1] = 1
  planedims[2] = 1
  if plane == 'xy':
    planedims[2] = 0
  elif plane == 'xz':
    planedims[1] = 0
  elif plane == 'yz':
    planedims[0] = 0
  elif plane != '3d':
    raise Exception('Invalid argument for specifying plane, must be xy, xz, or yz plane')

  for i in range(3):
    c_pos[i] = pos[i]

  il = <IntList*> malloc(sizeof(IntList))
  c_analyze.nbhood(c_pos, r_catch, il, planedims)

  result = create_nparray_from_IntList(il)
  realloc_intlist(il, 0)
  free(il)
  return result

#
# Pressure analysis
#
def pressure(self, pressure_type = 'all', id1 = 'default', id2 = 'default', v_comp=False):
  cdef vector[string] pressure_labels
  cdef vector[double] pressures

  checkTypeOrExcept(v_comp, 1, bool, "v_comp must be a boolean")

  if pressure_type=='all':
    c_analyze.analyze_pressure_all(pressure_labels, pressures, v_comp)
    return pressure_labels, pressures
  elif id1 == 'default' and id2 == 'default':
    pressure = c_analyze.analyze_pressure(pressure_type, v_comp)
    return pressure
  elif id1 != 'default' and id2 == 'default':
    checkTypeOrExcept(id1, 1, int, "id1 must be an int")
    pressure = c_analyze.analyze_pressure_single(pressure_type, id1, v_comp)
    return pressure
  else:
    checkTypeOrExcept(id1, 1, int, "id1 must be an int")
    checkTypeOrExcept(id2, 1, int, "id2 must be an int")
    pressure = c_analyze.analyze_pressure_pair(pressure_type, id1, id2, v_comp)
    return pressure

def stress_tensor(self, stress_type = 'all', id1 = 'default', id2 = 'default', v_comp=False):
  cdef vector[string] stress_labels
  cdef vector[double] stresses

  checkTypeOrExcept(v_comp, 1, bool, "v_comp must be a boolean")
  
  if stress_type=='all':
    c_analyze.analyze_stress_tensor_all(stress_labels, stresses, v_comp)
    return stress_labels, stresses
  elif id1 == 'default' and id2 == 'default':
    if (c_analyze.analyze_stress_tensor(stress_type, v_comp, stresses)):
      raise Exception("Error while calculating stress tensor")
    return stresses
  elif id1 != 'default' and id2 == 'default':
    checkTypeOrExcept(id1, 1, int, "id1 must be an int")
    if (c_analyze.analyze_stress_single(stress_type, id1, v_comp, stresses)):
      raise Exception("Error while calculating stress tensor")
    return stresses
  else:
    checkTypeOrExcept(id1, 1, int, "id1 must be an int")
    checkTypeOrExcept(id2, 1, int, "id2 must be an int")
    if (c_analyze.analyze_stress_pair(stress_type, id1, id2, v_comp, stresses)):
      raise Exception("Error while calculating stress tensor")
    return stresses
    
def local_stress_tensor(self, periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1)):
  cdef DoubleList* local_stress_tensor=NULL
  cdef int[3] c_periodicity, c_bins
  cdef double[3] c_range_start, c_stress_range

  for i in range(3):
    c_bins[i]=bins[i]
    c_periodicity[i]=periodicity[i]
    c_range_start[i]=range_start[i]
    c_stress_range[i]=stress_range[i]
  
  if c_analyze.analyze_local_stress_tensor(c_periodicity, c_range_start, c_stress_range, c_bins, local_stress_tensor):
    raise Exception("Error while calculating local stress tensor")
  stress_tensor =  create_nparray_from_DoubleList(local_stress_tensor)
  free (local_stress_tensor)
  return stress_tensor

#
# Energy analysis
#
def energy(system, etype = 'all', id1 = 'default', id2 = 'default'):
  if system.n_part == 0:
    print  'no particles'
    return 'no particles'

  if c_analyze.total_energy.init_status == 0:
    c_analyze.init_energies(&c_analyze.total_energy)
    c_analyze.master_energy_calc()
  _value = 0.0

  if etype == 'all':
    _result = energy(system, 'total') + ' ' + energy(system, 'kinetic')
    _result += energy(system, 'nonbonded',0,0)
    # todo: check for existing particle and bond types
    # and add those to _result
    return _result

  if etype == 'total':
    if id1 != 'default' or id2 != 'default':
      print ('warning: energy(\'total\') does not need '
             'further arguments, ignored.')
    for i in range(c_analyze.total_energy.data.n):
      _value += c_analyze.total_energy.data.e[i]
    return '{ energy: %f }' % _value

  if etype == 'kinetic':
    if id1 != 'default' or id2 != 'default':
      print ('warning: energy(\'kinetic\') does not need '
             'further arguments, ignored.')
    _value = c_analyze.total_energy.data.e[0]
    return '{ kinetic: %f }' % _value

  # coulomb interaction
  if etype == 'coulomb':
    if(code_info.electrostatics_defined()):
      for i in range(c_analyze.total_energy.n_coulomb):
        _value += c_analyze.total_energy.coulomb[i]
      return '{ coulomb: %f }' % _value
    else:
      print  'error: ELECTROSTATICS not compiled'
      return 'error: ELECTROSTATICS not compiled'

  if etype == 'magnetic':
    if(code_info.dipoles_defined()):
      for i in range(c_analyze.total_energy.n_dipolar):
        _value += c_analyze.total_energy.dipolar[i]
      return '{ magnetic: %f }' % _value
    else:
      print  'error: DIPOLES not compiled'
      return 'error: DIPOLES not compiled'

  # bonded interactions
  if etype == 'bonded':
    if not isinstance(id1, int):
      print ('error: analyze.energy(\'bonded\',<bondid>): '
             '<bondid> must be integer')
      raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                      '<bondid> must be integer')
    else:
    # todo: check if bond type id1 exist
      _value = c_analyze.obsstat_bonded(&c_analyze.total_energy, id1)[0]
      return '{ %d bonded: %f }' % (id1,_value)

  # nonbonded interactions
  if etype == 'nonbonded':
    if not isinstance(id1, int):
      print  ('error: analyze.energy(\'bonded\',<bondid>): '
              '<bondid> must be integer')
      raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                      '<bondid> must be integer')
    if not isinstance(id2, int):
      print  ('error: analyze.energy(\'bonded\',<bondid>): '
              '<bondid> must be integer')
      raise TypeError('analyze.energy(\'bonded\',<bondid>): '
                      '<bondid> must be integer')
    else:
    # todo: check if particle types id1 and id2 exist
      _value = c_analyze.obsstat_nonbonded(&c_analyze.total_energy, id1, id2)[0]
      return '{ %d %d nonbonded: %f }' % (id1,id2,_value)

  return 'error: unknown feature of analyze energy: \'%s\'' % etype
