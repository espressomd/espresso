import numpy
import os
if not "ETS_TOOLKIT" in os.environ:
	os.environ["ETS_TOOLKIT"] = "wx"
from mayavi import mlab
from pyface.api import GUI
from pyface.timer.api import Timer
from tvtk.tools import visual
import atexit
import threading

from espressomd.interactions import NonBondedInteractions

# workaround for https://github.com/enthought/mayavi/issues/3
import vtk
output=vtk.vtkFileOutputWindow()
output.SetFileName("/dev/null")
vtk.vtkOutputWindow().SetInstance(output)

class mayavi_live:
	"""This class provides live visualization using Enthought Mayavi.
	Use the update method to push your current simulation state after
	integrating. If you run your integrate loop in a separate thread, 
	you can call run_gui_event_loop in your main thread to be able to
	interact with the GUI."""
	def __init__(self, system):
		"""Constructor.
		**Arguments**

		:system: instance of espressomd.System
		"""
		self.system = system

		# objects drawn
		self.points = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="sphere", scale_factor=1, name="Particles")
		self.points.glyph.color_mode = 'color_by_scalar'
		self.points.glyph.glyph_source.glyph_source.center = [0, 0, 0]
		self.box = mlab.outline(extent=(0,0,0,0,0,0), color=(1,1,1), name="Box")
		self.arrows = mlab.quiver3d([],[],[], [],[],[], scalars=[], mode="2ddash", scale_factor=1, name="Bonds")
		self.arrows.glyph.color_mode = 'color_by_scalar'

		# state
		self.data = None
		self.last_N = 1
		self.last_Nbonds = 1
		self.last_boxl = [0,0,0]
		self.running = False
		self.last_T = None

		# GUI window
		self.gui = GUI()
		self.timers = [Timer(100, self._draw)]

	def _draw(self):
		"""Update the Mayavi objects with new particle information.
		This is called periodically in the GUI thread"""
		if self.data is None:
			return

		f = mlab.gcf()
		visual.set_viewer(f)

		coords, types, radii, N_changed, bonds, Nbonds_changed, boxl, box_changed = self.data
		self.data = None

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
		"""Pull the latest particle information from Espresso.
		This is the only function that should be called from the computation thread.
		It does not call any Mayavi functions unless it is being called from the main (GUI) thread."""

		if self.last_T is not None and self.last_T == self.system.time:
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

		if self.data is None:
			self.data = coords, types, radii, (self.last_N != N), \
			            bond_coords, (self.last_Nbonds != Nbonds), \
			            boxl, (self.last_boxl != boxl).any()
		else:
			self.data = coords, types, radii, self.data[3] or (self.last_N != N), \
			            bond_coords, self.data[5] or (self.last_Nbonds != Nbonds), \
			            boxl, self.data[7] or (self.last_boxl != boxl).any()
		self.last_N = N
		self.last_Nbonds = Nbonds
		self.last_boxl = boxl

		# when drawing from the main thread, the timer never fires, but we can safely call draw ourselves
		if isinstance(threading.current_thread(), threading._MainThread):
			self._draw()

	def process_gui_events(self):
		"""Process GUI events, e.g. mouse clicks, in the Mayavi window.
		Call this function as often as you can to get a smooth GUI experience."""
		self.gui.process_events()

	def run_gui_event_loop(self):
		"""Start the GUI event loop.
		This function blocks until the Mayavi window is closed.
		So you should only use it if your Espresso simulation's integrate loop is running in a secondary thread."""
		self.gui.start_event_loop()

	def register_callback(self, cb, interval=1000):
		self.timers.append(Timer(interval, cb))

# TODO: constraints
