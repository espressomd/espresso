import re
import numpy
import sys
import os
try:
    import __builtin__ # for Python 2
except ImportError:
    import builtins as __builtin__ # for Python 3

col_types = {'^di': 3,
 '^ext_f': 3,
 '^ext_t': 3,
 '^f': 3,
 '^fix': 3,
 '^i': 1,
 '^ma': 1,
 '^mo': 1,
 '^omega': 3,
 '^omega_b': 3,
 '^omega_l': 3,
 '^p': 3,
 '^q$': 1,
 '^qu': 4,
 '^tbf': 3,
 '^torque': 3,
 '^torque_b': 3,
 '^torque_l': 3,
 '^ty': 1,
 '^v': 3,
 '^vi': 1,
 '^vs': 2}

re_block_type = re.compile('^{([a-z_]+)\s')
re_particles = re.compile('{particles\s*{([a-z_\s]*)}((?:\s*{.*?})*)\s*}')
re_particle = re.compile('\{(\s*(?:[0-9\.-]+\s*)*)\}')
re_space = re.compile('\s')
re_variable = re.compile('\{(?!variable\s)([^\s\}]+)([^\}]*)\}')
re_int_list = re.compile('^[0-9\s]+$')
re_float_list = re.compile('^[0-9\.\s]+$')
re_int = re.compile('[0-9]+')
re_float = re.compile('[0-9\.]+')

re_col_types = {}
for col_type in col_types:
	re_col_types[col_type] = re.compile(col_type)

def load_col_types(blockfile_support_tcl):
	global col_types, re_col_types
	col_types, re_col_types = {}, {}
	with __builtin__.open(blockfile_support_tcl) as f:
		for i in re.finditer('"(\^[a-z_]+\$*)"\s*{.*?; incr idx ([0-9]*)\s*}', f.read(), re.DOTALL):
			var_name, var_cols = i.group(1), i.group(2)
			if var_cols == '': var_cols = 1
			else: var_cols = int(var_cols)
			col_types[var_name] = var_cols
	for col_type in col_types:
		re_col_types[col_type] = re.compile(col_type)

def process(block):
	"""
	Processes the block and returns the block's type and the block's contents as a tuple.
	"""
	block_type = re_block_type.match(block).group(1)
	
	# particle data looks like this:
	# {particles {id pos v q f} 
	# 	{0 7.875 9.0 0.0 0.0 0.0 0.0 -0.6 7.29706 -1.08036 -3.39398}
	# 	{1 7.875 9.0 1.2 0.0 0.0 0.0 -0.6 1.37885 -0.320172 -2.49209}
	if block_type == 'particles':
		particles = re_particles.match(block)
		
		# get the field labels and determine how many columns a field consists of (i.e. dimensionality of vector quantities)
		fields = particles.group(1).split(' ')
		field_lengths = [0]*len(fields)
		for i,field in enumerate(fields):
			for col_type in col_types:
				if re_col_types[col_type].match(field) is not None:
					field_lengths[i] = col_types[col_type]
					break
			if field_lengths[i] == 0:
				raise Exception("Field '%s' is unknown. Dimensionality of this field cannot be determined." % field)
		
		# split the particle block into the individual particles
		particles = [re_space.split(x) for x in re_particle.findall(particles.group(2))]
		N = len(particles)
		
		# create one empty numpy array per field with the appropriate dimensionality
		particle_data = {}
		for i,field in enumerate(fields):
			dtype = float
			if re_col_types['^i'].match(field) is not None or re_col_types['^ty'].match(field) is not None:
				# the ID and type fields are integer, everything else is a float
				dtype = int
			particle_data[field] = numpy.empty((N,field_lengths[i]), dtype=dtype)
		
		# fill the numpy arrays with the particle data
		for i,particle in enumerate(particles):
			fieldcount = 0
			for j,field in enumerate(fields):
				particle_data[field][i] = particle[fieldcount:fieldcount+field_lengths[j]]
				fieldcount += field_lengths[j]
		
		return block_type, particle_data
	
	# variables look like this:
	# {variable  {box_l 18.0 18.0 480.0} {test 18.0} {stringtest test} }
	if block_type == 'variable':
		variables = {}
		# extract all variable names and values from the block
		for m in re_variable.finditer(block):
			try:
				# convert it to an integer/float if it's a (list of) integer/float variable(s)
				value = m.group(2).strip()
				if re_int_list.match(value) is not None:
					value = numpy.array(re_int.findall(value), dtype=int)
				elif re_float_list.match(value) is not None:
					value = numpy.array(re_float.findall(value), dtype=float)
				variables[m.group(1)] = value
			except:
				# return it as a string if we can't convert it to a number
				variables[m.group(1)] = m.group(2).strip()
		
		return block_type, variables
	
	# any other block types are returned as string
	return block_type, block[len(block_type)+1:-2].strip()

class blockfile(object):
	f = None
	
	def __init__(self, path):
		"""
		Opens the blockfile.
		"""
		self.f = __builtin__.open(path)
	
	def __iter__(self):
		"""
		Iterator over all the blocks in the open blockfile.
		The iterator returns each block as a tuple of block type and block contents.
		"""
		block = ''
		while True:
			block += self.f.readline() # this assumes that all blocks end with a newline
			if block.count('{') == block.count('}'): # we have a full block:
				yield process(block)
				block = ''
			if self.f.tell() == os.fstat(self.f.fileno()).st_size:
				return
	
	def __del__(self):
		"""
		Closes the blockfile.
		"""
		self.close()
	
	def close(self):
		"""
		Closes the blockfile.
		"""
		self.f.close()

if __name__ != "__main__":
	def open(path):
		"""
		Opens the blockfile at path.
		"""
		return blockfile(path)

if __name__ == "__main__":
	# If this script is called directly with the argument --extract, it extracts the list of column types from scripts/blockfile_support.tcl.
	if len(sys.argv) > 1 and sys.argv[1] == '--extract':
		if len(sys.argv) < 2 or not os.path.exists(sys.argv[2]):
			sys.stderr.write("Please specify the path to scripts/blockfile_support.tcl from Espresso.\n")
			sys.exit(1)
		load_col_types(sys.argv[2])
		from pprint import pprint
		print(pprint(col_types))
		sys.exit()
