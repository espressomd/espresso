import numpy
cimport numpy
import os
from libcpp cimport bool
from espressomd.particle_data import ParticleHandle
from particle_data cimport *
from espressomd.interactions cimport *
from espressomd.system cimport *
from libcpp.vector cimport vector

include "myconfig.pxi"

cdef extern from "utils.hpp":
    void get_mi_vector(double * res, double * a, double * b)

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

cdef class mayaviLive(object):
    """
    This class provides live visualization using Enthought Mayavi.  Use the
    update method to push your current simulation state after integrating. If
    you run your integrate loop in a separate thread, you can call
    run_gui_event_loop in your main thread to be able to interact with the GUI.
    
    Parameters
    ----------
    system : instance of espressomd.System
    particle_sizes : (optional) function, list, or dict, which maps particle types to radii

    """

    cdef object system
    cdef object particle_sizes
    cdef object points 
    cdef object box 
    cdef object arrows
    cdef object data 
    cdef int last_N 
    cdef int last_Nbonds 
    cdef double last_boxl[3] 
    cdef bool running 
    cdef double last_T 

    cdef object gui 
    cdef object timers 


    def __init__(self, system, particle_sizes='auto'):
        self.system = system
        self.particle_sizes = particle_sizes

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
        self.last_Nbonds = 0
        self.last_boxl = [0,0,0]
        self.running = False
        self.last_T = -1

        # GUI window
        self.gui = GUI()
        self.timers = [Timer(100, self._draw)]

    def _determine_radius(self, t):
        """Determine radius of particle type t for visualization."""
        def radius_from_lj(t):
            radius = 0.
            IF LENNARD_JONES:
                try:
                    radius = 0.5 * get_ia_param(t,t).LJ_sig
                except:
                    radius = 0.

            if radius == 0:
                radius = 0.5  # fallback value
            return radius

        radius = 0.
        if self.particle_sizes is 'auto':
            radius = radius_from_lj(t)
        elif callable(self.particle_sizes):
            radius = self.particle_sizes(t)
        else:
            try:
                radius = self.particle_sizes[t]
            except:
                radius = radius_from_lj(t)
        return radius

    def _draw(self):
        """
        Update the Mayavi objects with new particle information.
        This is called periodically in the GUI thread.
        
        """
        if self.data is None:
            return

        assert isinstance(threading.current_thread(), threading._MainThread)

        f = mlab.gcf()
        visual.set_viewer(f)

        coords, types, radii, N_changed, bonds, Nbonds_changed, boxl, box_changed = self.data
        self.data = None

        if box_changed or not self.running:
            self.box.set(bounds=(0,boxl[0], 0,boxl[1], 0,boxl[2]))
        if not N_changed:
            self.points.mlab_source.set(x=coords[:,0]%boxl[0], y=coords[:,1]%boxl[1], z=coords[:,2]%boxl[2], u=radii, v=radii, w=radii, scalars=types)
        else:
            self.points.mlab_source.reset(x=coords[:,0]%boxl[0], y=coords[:,1]%boxl[1], z=coords[:,2]%boxl[2], u=radii, v=radii, w=radii, scalars=types)
        if not self.running:
            f.scene.reset_zoom()
            self.running = True

        if not Nbonds_changed:
            if bonds.shape[0] > 0:
                self.arrows.mlab_source.set(x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5])
        else:
            self.arrows.mlab_source.reset(x=bonds[:,0], y=bonds[:,1], z=bonds[:,2], u=bonds[:,3], v=bonds[:,4], w=bonds[:,5], scalars=bonds[:,6])

    def update(self):
        """Pull the latest particle information from Espresso.
        
        This is the only function that should be called from the computation
        thread.  It does not call any Mayavi functions unless it is being
        called from the main (GUI) thread.
        
        """

        if self.last_T is not None and self.last_T == self.system.time:
            return
        self.last_T = self.system.time

        cdef int N = len(self.system.part)
        coords = numpy.zeros((N,3))
        types = numpy.empty(N, dtype=int)
        inter = NonBondedInteractions()
        radii = numpy.empty(N)
        cdef int i=0,j=0,k=0 
        cdef int t 
        cdef int partner
        cdef const particle* p
        cdef IA_parameters* ia
        cdef vector[int] bonds

    # Using (additional) untyped variables and python constructs in the loop 
    # will slow it down considerably. 
        for i in range(N):
            p = get_particle_data_ptr(i)
            if not p:
                continue

            coords[j,:] = numpy.array([p.r.p[0],p.r.p[1],p.r.p[2]])
            t = p.p.type
            types[j] = t +1
            radii[j] = self._determine_radius(t)

            # Iterate over bonds
            k = 0
            while k<p.bl.n:
                # Bond type
                t = p.bl.e[k]
                k += 1
                # Iterate over bond partners and store each connection
                for l in range(bonded_ia_params[t].num):
                    bonds.push_back(i)
                    bonds.push_back(p.bl.e[k])
                    bonds.push_back(t)
                    k += 1
            j += 1
        assert j == len(self.system.part)
        cdef int Nbonds = bonds.size()//3

        bond_coords = numpy.empty((Nbonds,7))

        cdef int n
        cdef const particle* p1
        cdef const particle* p2
        cdef double bond_vec[3]
        for n in range(Nbonds):
            i = bonds[3*n]
            j = bonds[3*n+1]
            t = bonds[3*n+2]
            p1 = get_particle_data_ptr(i)
            p2 = get_particle_data_ptr(j)
            bond_coords[n,:3] = numpy.array([p1.r.p[0],p1.r.p[1],p1.r.p[2]])
            get_mi_vector(bond_vec, p2.r.p.data(), p1.r.p.data())
            bond_coords[n,3:6] = bond_vec
            bond_coords[n,6] = t

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

    def processGuiEvents(self):
        """Process GUI events, e.g. mouse clicks, in the Mayavi window.
        
        Call this function as often as you can to get a smooth GUI experience.
        
        """
        assert isinstance(threading.current_thread(), threading._MainThread)
        self.gui.process_events()

    def start(self):
        """Start the GUI event loop.
        
        This function blocks until the Mayavi window is closed. So you should
        only use it if your Espresso simulation's integrate loop is running in
        a secondary thread.
        
        """
        assert isinstance(threading.current_thread(), threading._MainThread)
        self.gui.start_event_loop()

    def register_callback(self, cb, interval=1000):
        self.timers.append(Timer(interval, cb))

# TODO: constraints
