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
include "myconfig.pxi"
cimport c_analyze
cimport utils
cimport particle_data
import utils
import code_info
import particle_data
from libcpp.string cimport string #import std::string as string
from libcpp.vector cimport vector #import std::vector as vector
from interactions import *
from interactions cimport *

#
# Minimal distance between particles
#

def mindist(system, p1 = 'default', p2 = 'default'):
  """Minimal distance between particles
    mindist(p1="default",p2="default")
    
    p1, p2: lists of particle types
  """

  cdef IntList* set1
  cdef IntList* set2

  if p1 == 'default' and p2 == 'default':
    result = c_analyze.mindist(NULL,NULL)
  elif (p1 == 'default' and not p2 == 'default') or \
       (not p1 == 'default' and p2 == 'default'):
    raise Exception("Both, p1 and p2 have to be specified\n"+__doc__)
  else:
    for i in range(len(p1)):
      if not isinstance(p1[i],int):
        raise ValueError("Particle types in p1 and p2 have to be of type int: "+str(p1[i]))

    for i in range(len(p2)):
      if not isinstance(p2[i],int):
        raise ValueError("Particle types in p1 and p2 have to be of type int"+str(p2[i]))
  
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
def distto(system, id=None,pos=None):
  """"Distance to particle or point
      distto(id=None,pos=None)
      id: id of the particle
      pos: position as array
  """

  if id==None and pos==None:
    raise Exception("Either id or pos have to be specified\n"+__doc__)
  
  if id !=None and pos != None:
    raise Exception("Only one of id or pos may be specified\n"+__doc__)


  cdef double cpos[3]
  if system.n_part == 0:
    raise Exception("no particles")

  
  # Get position 
  # If particle id specified
  if id !=None:
    if not isinstance(id,int):
      raise ValueError("Id has to be an integer")
    _pos = system.part[id].pos
    for i in range(3):
      cpos[i] = _pos[i]
  else:
    for i in range(3):
      cpos[i] = pos[i]
      _id = -1
  return c_analyze.distto(cpos,_id)

# get all particles in neighborhood r_catch of pos and return their ids
# in il. plane can be used to specify the distance in the xy, xz or yz
# plane
def nbhood(self, pos=None, r_catch=None, plane = '3d'):
  """nbhood(pos=None, r_catch=None, plane = '3d'):
     get all particles in neighborhood r_catch of pos and return their ids
     in il. plane can be used to specify the distance in the xy, xz or yz
     plane
     nbhood(pos=None, r_catch=None, plane = '3d'):
  """
     
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
def pressure(self,v_comp=0):
  """Pewaauew
     pressure(pressure_type = 'all', id1 = 'default', id2 = 'default', v_comp=False)
  """
#
  checkTypeOrExcept(v_comp, 1, int, "v_comp must be a boolean")
#
  # Dict to store the results
  p={}

  # Update in espresso core if necessary
  if (c_analyze.total_pressure.init_status != 1+v_comp):
    c_analyze.update_pressure(v_comp)
#
  # Individual components of the pressure

  # Total pressure
  cdef int i
  cdef double tmp
  tmp=0
  for i in range(c_analyze.total_pressure.data.n):
    tmp+=c_analyze.total_pressure.data.e[i]

  p["total"]=tmp

  # Ideal
  p["ideal"] =c_analyze.total_pressure.data.e[0]

  # Nonbonded
  cdef double total_bonded
  total_bonded=0
  for i in range(c_analyze.n_bonded_ia):
    if (bonded_ia_params[i].type != 0): 
     p["bonded",i] = c_analyze.obsstat_bonded(&c_analyze.total_pressure, i)[0]
     total_bonded += c_analyze.obsstat_bonded(&c_analyze.total_pressure, i)[0]
  p["bonded"] = total_bonded

  # Non-Bonded interactions, total as well as intra and inter molecular
  cdef int j
  cdef double total_intra
  cdef double total_inter
  cdef double total_non_bonded
  total_inter=0
  total_intra=0
  total_non_bonded=0

  for i in range(c_analyze.n_particle_types):
    for j in range(c_analyze.n_particle_types):
#      if checkIfParticlesInteract(i, j):
        p["nonBonded",i,j] =c_analyze.obsstat_nonbonded(&c_analyze.total_pressure, i, j)[0]
        total_non_bonded =c_analyze.obsstat_nonbonded(&c_analyze.total_pressure, i, j)[0]
        total_intra +=c_analyze.obsstat_nonbonded_intra(&c_analyze.total_pressure_non_bonded, i, j)[0]
        p["nonBondedIntra",i,j] =c_analyze.obsstat_nonbonded_intra(&c_analyze.total_pressure_non_bonded, i, j)[0]
        p["nonBondedInter",i,j] =c_analyze.obsstat_nonbonded_inter(&c_analyze.total_pressure_non_bonded, i, j)[0]
        total_inter+= c_analyze.obsstat_nonbonded_inter(&c_analyze.total_pressure_non_bonded, i, j)[0]
  p["nonBondedIntra"]=total_intra
  p["nonBondedInter"]=total_inter
  p["nonBondedInter"]=total_inter
  p["nonBonded"]=total_non_bonded

  # Electrostatics
  IF ELECTROSTATICS == 1:
    cdef double total_coulomb
    total_coulomb=0
    for i in range(c_analyze.total_pressure.n_coulomb):
      total_coulomb+=c_analyze.total_pressure.coulomb[i]
      p["coulomb",i]=c_analyze.total_pressure.coulomb[i]
    p["coulomb"]=total_coulomb

  # Dipoles
  IF DIPOLES == 1:
    cdef double total_dipolar
    total_dipolar=0
    for i in range(c_analyze.total_pressure.n_dipolar):
      total_dipolar+=c_analyze.total_pressure.dipolar[i]
      p["dipolar",i]=c_analyze.total_pressure.coulomb[i]
    p["dipolar"]=total_dipolar

  # virtual sites
  IF VIRTUAL_SITES_RELATIVE ==1:
    p["vs_relative"]=c_analyze.total_pressure.vs_relative[0]

  return p


def stress_tensor(self, stress_type = 'all', id1 = 'default', id2 = 'default', v_comp=False):
  """stress_tensor(self, stress_type = 'all', id1 = 'default', id2 = 'default', v_comp=False)"""
  cdef vector[string] stress_labels
  cdef vector[double] stresses

  checkTypeOrExcept(v_comp, 1, bool, "v_comp must be a boolean")
  
  if stress_type=='all':
    c_analyze.analyze_stress_tensor_all(stress_labels, stresses, v_comp)
    return dict(zip(stress_labels, stresses))
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
  """local_stress_tensor(periodicity=(1, 1, 1), range_start=(0.0, 0.0, 0.0), stress_range=(1.0, 1.0, 1.0), bins=(1, 1, 1))
  """

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
  """energy(system, etype = 'all', id1 = 'default', id2 = 'default')"""
#  if system.n_part == 0:
#    raise Exception('no particles')

  e={}

  if c_analyze.total_energy.init_status == 0:
    c_analyze.init_energies(&c_analyze.total_energy)
    c_analyze.master_energy_calc()
  
 
  # Individual components of the pressur

  # Total energy 
  cdef int i
  cdef double tmp
  tmp=0
  for i in range(c_analyze.total_energy.data.n):
    tmp+=c_analyze.total_energy.data.e[i]

  e["total"]=tmp

  # Ideal
  e["ideal"] =c_analyze.total_energy.data.e[0]

  # Nonbonded
  cdef double total_bonded
  total_bonded=0
  for i in range(c_analyze.n_bonded_ia):
    if (bonded_ia_params[i].type != 0): 
     e["bonded",i] = c_analyze.obsstat_bonded(&c_analyze.total_energy, i)[0]
     total_bonded += c_analyze.obsstat_bonded(&c_analyze.total_energy, i)[0]
  e["bonded"] = total_bonded

  # Non-Bonded interactions, total as well as intra and inter molecular
  cdef int j
  cdef double total_intra
  cdef double total_inter
  cdef double total_non_bonded
  total_inter=0
  total_intra=0
  total_non_bonded=0

  for i in range(c_analyze.n_particle_types):
    for j in range(c_analyze.n_particle_types):
#      if checkIfParticlesInteract(i, j):
        e["nonBonded",i,j] =c_analyze.obsstat_nonbonded(&c_analyze.total_energy, i, j)[0]
        total_non_bonded =c_analyze.obsstat_nonbonded(&c_analyze.total_energy, i, j)[0]
#        total_intra +=c_analyze.obsstat_nonbonded_intra(&c_analyze.total_energy_non_bonded, i, j)[0]
#        e["nonBondedIntra",i,j] =c_analyze.obsstat_nonbonded_intra(&c_analyze.total_energy_non_bonded, i, j)[0]
#        e["nonBondedInter",i,j] =c_analyze.obsstat_nonbonded_inter(&c_analyze.total_energy_non_bonded, i, j)[0]
#        total_inter+= c_analyze.obsstat_nonbonded_inter(&c_analyze.total_energy_non_bonded, i, j)[0]
#  e["nonBondedIntra"]=total_intra
#  e["nonBondedInter"]=total_inter
  e["nonBonded"]=total_non_bonded

  # Electrostatics
  IF ELECTROSTATICS == 1:
    cdef double total_coulomb
    total_coulomb=0
    for i in range(c_analyze.total_energy.n_coulomb):
      total_coulomb+=c_analyze.total_energy.coulomb[i]
      e["coulomb",i]=c_analyze.total_energy.coulomb[i]
    e["coulomb"]=total_coulomb

  # Dipoles
  IF DIPOLES == 1:
    cdef double total_dipolar
    total_dipolar=0
    for i in range(c_analyze.total_energy.n_dipolar):
      total_dipolar+=c_analyze.total_energy.dipolar[i]
      e["dipolar",i]=c_analyze.total_energy.coulomb[i]
    e["dipolar"]=total_dipolar

  return e 



