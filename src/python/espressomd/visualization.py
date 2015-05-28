import numpy
from mayavi import mlab
from multiprocessing import Process, Queue
import atexit

from espressomd.interactions import NonBondedInteractions

# workaround for https://github.com/enthought/mayavi/issues/3
import vtk
output=vtk.vtkFileOutputWindow()
output.SetFileName("/dev/null")
vtk.vtkOutputWindow().SetInstance(output)

def mayavi_render(queue):
	points = mlab.points3d([0],[0],[0])
	arrows = mlab.quiver3d([0],[0],[0], [0],[0],[0])

	@mlab.show
	@mlab.animate(delay=10, ui=True)
	def animate():
		f = mlab.gcf()
		running = False
		while True:
			coords, types, radii, N_changed, bonds, Nbonds_changed, boxl = queue.get()

			if not N_changed:
				points.mlab_source.set(x=coords[:,0], y=coords[:,1], z=coords[:,2])
			else:
				points.mlab_source.reset(x=coords[:,0], y=coords[:,1], z=coords[:,2])
			if not running:
				f.scene.reset_zoom()
				running = True

			if not Nbonds_changed:
				arrows.mlab_source.set  (x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5])
			else:
				arrows.mlab_source.reset(x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5])

			print boxl

			yield
	animate()

class mayavi_live:
	def __init__(self, system):
		self.system = system
		self.queue = Queue()
		self.process = Process(target=mayavi_render, args=(self.queue,))
		self.process.start()
		self.last_N = 1
		self.last_Nbonds = 1
		atexit.register(self.process.terminate)

	def __del__(self):
		self.process.terminate()

	def update(self):
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
			radii[i] = inter[t,t].lennardJones._getParamsFromEsCore()['sigma'] * 0.5

			bs = self.system.part[i].bonds
			for b in bs:
				for p in b[1:]:
					bonds.append((i,p))
		Nbonds = len(bonds)
		bond_coords = numpy.empty((Nbonds+1,6))
		bond_coords[0] = [0,0,0,0,0,0] # don't pass an empty array to Mayavi

		for n in range(Nbonds):
			i,j = bonds[n]
			bond_coords[n+1,:3] = self.system.part[i].pos
			bond_coords[n+1,3:] = self.system.part[j].pos - self.system.part[i].pos

		boxl = self.system.box_l
		print boxl, type(boxl)

		self.queue.put((coords, types, radii, (self.last_N != N), 
		                bond_coords, (self.last_Nbonds != Nbonds),
		                boxl))
		self.last_N = N
		self.last_Nbonds = Nbonds

# TODO: show simulation box: http://docs.enthought.com/mayavi/mayavi/auto/example_mlab_visual.html
# TODO: color particles by id and set their radius to sigma/2: http://docs.enthought.com/mayavi/mayavi/mlab_changing_object_looks.html
# TODO: constraints
# CHECK: bonds
