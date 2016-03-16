import sys
if sys.platform == 'darwin':
	try:
		import multiprocessing
		multiprocessing.set_start_method('spawn')
	except:
		raise Exception("Mayavi visualization is not supported on Mac OS X because fork()ed processes may not have a GUI.")

import numpy
from mayavi import mlab
from tvtk.tools import visual
from multiprocessing import Process, Queue
import atexit

from espressomd.interactions import NonBondedInteractions

# workaround for https://github.com/enthought/mayavi/issues/3
import vtk
output=vtk.vtkFileOutputWindow()
output.SetFileName("/dev/null")
vtk.vtkOutputWindow().SetInstance(output)

def mayavi_render(queue):
	points = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="sphere", scale_factor=1, name="Particles")
	points.glyph.color_mode = 'color_by_scalar'
	points.glyph.glyph_source.glyph_source.center = [0, 0, 0]
	box = mlab.outline(extent=(0,0,0,0,0,0), color=(1,1,1), name="Box")
	arrows = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="2ddash", scale_factor=1, name="Bonds")
	arrows.glyph.color_mode = 'color_by_scalar'

	def animate():
		f = mlab.gcf()
		visual.set_viewer(f)

		running = False
		while True:
			coords, types, radii, N_changed, bonds, Nbonds_changed, boxl, box_changed = queue.get()

			if box_changed or not running:
				box.set(bounds=(0,boxl[0], 0,boxl[1], 0,boxl[2]))

			if not N_changed:
				points.mlab_source.set(x=coords[:,0], y=coords[:,1], z=coords[:,2], u=radii, v=radii, w=radii, scalars=types)
			else:
				points.mlab_source.reset(x=coords[:,0], y=coords[:,1], z=coords[:,2], u=radii, v=radii, w=radii, scalars=types)
			if not running:
				f.scene.reset_zoom()
				running = True

			if not Nbonds_changed:
				if bonds.shape[0] > 0:
					arrows.mlab_source.set  (x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5])
			else:
				arrows.mlab_source.reset(x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5], scalars=bonds[:,6])

			yield

	# Remove the parameter window
	from mayavi.tools.animator import Animator
	def e(self):
		raise Exception("Please ignore this exception. It is needed to keep the Parameters window from showing up.")
	Animator.show = e

	mlab.show(mlab.animate(delay=10, ui=True)(animate))()
	raise Exception("Animation window closed")

class mayavi_live:
	def __init__(self, system):
		self.system = system
		self.queue = Queue()
		self.process = Process(target=mayavi_render, args=(self.queue,))
		self.process.daemon = True
		self.process.start()
		self.last_N = 1
		self.last_Nbonds = 1
		self.last_boxl = [0,0,0]
		atexit.register(self.stop)

	def stop(self):
		if self.process.is_alive():
			self.process.terminate()

	def __del__(self):
		self.stop()

	def update(self):
		if not self.process.is_alive():
			raise Exception("Animation process terminated")
		N = self.system.n_part
		coords = numpy.empty((N,3))
		types = numpy.empty(N, dtype=int)
		inter = NonBondedInteractions()
		radii = numpy.empty(N)
		bonds = []

		for i in range(N):
			coords[i,:] = self.system.part[i].pos
			t = self.system.part[i].type
			types[i] = t +1
			radii[i] = inter[t,t].lennard_jones.get_params()['sigma'] * 0.5

			bs = self.system.part[i].bonds
			for b in bs:
				t = b[0]
				for p in b[1:]:
					bonds.append((i,p,t))
		Nbonds = len(bonds)
		bond_coords = numpy.empty((Nbonds,7))

		for n in range(Nbonds):
			i,j,t = bonds[n]
			bond_coords[n,:3] = self.system.part[i].pos
			bond_coords[n,3:6] = self.system.part[j].pos - self.system.part[i].pos
			bond_coords[n,6] = 0 # numeric bond IDs have been removed

		boxl = self.system.box_l

		self.queue.put(( coords, types, radii, (self.last_N != N), 
		                 bond_coords, (self.last_Nbonds != Nbonds),
		                 boxl, (self.last_boxl != boxl).any() ))
		self.last_N = N
		self.last_Nbonds = Nbonds
		self.last_boxl = boxl

# TODO: constraints
