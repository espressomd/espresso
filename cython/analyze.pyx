# For C-extern Analysis

cimport c_analyze
cimport utils
cimport particle_data
import utils
import code_info
import global_variables
import particle_data

#
# Minimal distance between particles
#
def mindist(p1 = 'default', p2 = 'default'):

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
def distto(id_or_pos):
  cdef double cpos[3]
  if global_variables.GlobalsHandle().n_part == 0:
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

#
# Energy analysis
#
def energy(etype = 'all', id1 = 'default', id2 = 'default'):

  if global_variables.GlobalsHandle().n_part == 0:
    print  'no particles'
    return 'no particles'

  if c_analyze.total_energy.init_status == 0:
    c_analyze.init_energies(&c_analyze.total_energy)
    c_analyze.master_energy_calc()
  _value = 0.0

  if etype == 'all':
    _result = energy('total') + ' ' + energy('kinetic')
    _result += energy('nonbonded',0,0)
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
