# For C-extern Analysis

cimport c_analyze
cimport utils
import utils

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

def energy(etype = 'total'):

	if c_analyze.total_energy.init_status == 0:
		c_analyze.init_energies(&c_analyze.total_energy)
	c_analyze.master_energy_calc()
	value = 0.0

	if etype == 'total':
		for i in range(c_analyze.total_energy.data.n):
			value += c_analyze.total_energy.data.e[i]
		return "Total Energy: %f" % value

	if etype == 'kinetic':
		value = c_analyze.total_energy.data.e[0]
		return "Kinetic Energy: %f" % value
