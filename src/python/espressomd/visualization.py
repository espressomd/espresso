import numpy
from mayavi import mlab
from pyface.api import GUI
from tvtk.tools import visual
import atexit

from espressomd.interactions import NonBondedInteractions

# workaround for https://github.com/enthought/mayavi/issues/3
import vtk
output=vtk.vtkFileOutputWindow()
output.SetFileName("/dev/null")
vtk.vtkOutputWindow().SetInstance(output)

class mayavi_live:
	def __init__(self, system):
		self.system = system

		# objects drawn
		self.points = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="sphere", scale_factor=1, name="Particles")
		self.points.glyph.color_mode = 'color_by_scalar'
		self.points.glyph.glyph_source.glyph_source.center = [0, 0, 0]
		self.box = mlab.outline(extent=(0,0,0,0,0,0), color=(1,1,1), name="Box")
		self.arrows = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="2ddash", scale_factor=1, name="Bonds")
		self.arrows.glyph.color_mode = 'color_by_scalar'

		# state
		self.last_N = 1
		self.last_Nbonds = 1
		self.last_boxl = [0,0,0]
		self.running = False
		self.last_T = -1

		# GUI window
		self.gui = GUI()
		self.gui.invoke_later(self.update)

	def draw(self, coords, types, radii, N_changed, bonds, Nbonds_changed, boxl, box_changed):
		f = mlab.gcf()
		visual.set_viewer(f)

		if box_changed or not self.running:
			self.box.set(bounds=(0,boxl[0], 0,boxl[1], 0,boxl[2]))

		if not N_changed:
			self.points.mlab_source.set(x=coords[:,0], y=coords[:,1], z=coords[:,2], u=radii, v=radii, w=radii, scalars=types)
		else:
			self.points.mlab_source.reset(x=coords[:,0], y=coords[:,1], z=coords[:,2], u=radii, v=radii, w=radii, scalars=types)
		if not self.running:
			f.scene.reset_zoom()
			self.running = True

		if not Nbonds_changed:
			if bonds.shape[0] > 0:
				self.arrows.mlab_source.set  (x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5])
		else:
			self.arrows.mlab_source.reset(x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5], scalars=bonds[:,6])

	def update(self):
		self.process_gui_events()

		if self.last_T == self.system.time:
			return
		self.last_T = self.system.time
		
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

		self.draw( coords, types, radii, (self.last_N != N), 
		           bond_coords, (self.last_Nbonds != Nbonds),
		           boxl, (self.last_boxl != boxl).any() )
		self.last_N = N
		self.last_Nbonds = Nbonds
		self.last_boxl = boxl

		self.gui.invoke_later(self.update)

	def process_gui_events(self):
		self.gui.process_events()

	def run_gui_event_loop(self):
		self.gui.start_event_loop()

# TODO: constraints
