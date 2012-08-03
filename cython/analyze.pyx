# For C-extern Analysis

cimport c_analyze
cimport utils
import utils
import code_info
import global_variables

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
#	free (set1)
#	free (set2)

	return result

def energy(etype = 'all', id1 = 'default', id2 = 'default'):

#	if global_variables.GlobalsHandle().n_part == 0:
#		print  'no particles'
#		return 'no particles'

	if c_analyze.total_energy.init_status == 0:
		c_analyze.init_energies(&c_analyze.total_energy)
		c_analyze.master_energy_calc()
	value = 0.0

	if etype == 'all':
		return energy('total') + ' ' + energy('kinetic')

	if etype == 'total':
		if id1 != 'default' or id2 != 'default':
			print 'warning: energy(\'total\') does not need further arguments, ignored.'
		for i in range(c_analyze.total_energy.data.n):
			value += c_analyze.total_energy.data.e[i]
		return '{ energy: %f }' % value

	if etype == 'kinetic':
		if id1 != 'default' or id2 != 'default':
			print 'warning: energy(\'kinetic\') does not need further arguments, ignored.'
		value = c_analyze.total_energy.data.e[0]
		return '{ kinetic: %f }' % value

# coulomb interaction
	if etype == 'coulomb':
		if(code_info.electrostatics_defined()):
			for i in range(c_analyze.total_energy.n_coulomb):
				value += c_analyze.total_energy.coulomb[i]
			return '{ coulomb: %f }' % value
		else:
			print  'error: ELECTROSTATICS not compiled'
			return 'error: ELECTROSTATICS not compiled'

	if etype == 'magnetic':
		if(code_info.dipoles_defined()):
			for i in range(c_analyze.total_energy.n_dipolar):
				value += c_analyze.total_energy.dipolar[i]
			return '{ magnetic: %f }' % value
		else:
			print  'error: DIPOLES not compiled'
			return 'error: DIPOLES not compiled'

# bonded interactions
	if etype == 'bonded':
		if not isinstance(id1, int):
			print  'error: analyze.energy(\'bonded\',<bondid>): <bondid> must be integer'
			raise TypeError('analyze.energy(\'bonded\',<bondid>): <bondid> must be integer')
		else:
#			value = c_analyze->obsstat_bonded(&total_energy, id1)
			return '{ bonded: %f }' % value

# bonded interactions
	if etype == 'nonbonded':
		print 'Not implemented yet'

	return 'error: unknown feature of analyze energy: \'%s\'' % etype





