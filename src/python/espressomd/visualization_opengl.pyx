from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from math import *
import math
import numpy as np
import os
import time
import espressomd
import collections
import sys
from threading import Thread

import scipy.spatial
include "myconfig.pxi"

from copy import deepcopy


class openGLLive(object):

    """This class provides live visualization using pyOpenGL.
    Use the update method to push your current simulation state after
    integrating. Modify the appearance with a list of keywords.
    Timed callbacks can be registered via the register_callback method.
    Keyboad callbacks via  keyboardManager.register_button().

    Attributes
    ----------

    system : instance of :attr:`espressomd.System`
    window_size : array_like :obj:`int`, optional
                  Size of the visualizer window in pixels.
    name : :obj:`str`, optional
           The name of the visualizer window.
    background_color : array_like :obj:`int`, optional
                       RGB of the background.
    periodic_images : array_like :obj:`int`, optional
                      Periodic repetitions on both sides of the box in xyzdirection.
    draw_box : :obj:`bool`, optional
               Draw wireframe boundaries.
    draw_axis : :obj:`bool`, optional
                Draws xyz system axes.
    quality_particles : :obj:`int`, optional
                        The number of subdivisions for particle spheres.
    quality_bonds : :obj:`int`, optional
                    The number of subdivisions for cylindrical bonds.
    quality_arrows : :obj:`int`, optional
                     The number of subdivisions for external force arrows.
    arrows_radius : :obj:`float`, optional
                    Arrow radius for velocity, force and ext_force arrows
    quality_constraints : :obj:`int`, optional
                          The number of subdivisions for primitive constraints.
    close_cut_distance : :obj:`float`, optional
                         The distance from the viewer to the near clipping plane.
    far_cut_distance : :obj:`float`, optional
                       The distance from the viewer to the far clipping plane.
    camera_position : :obj:`str` or array_like :obj:`float`, optional
                      Initial camera position. ``auto`` (default) for shiftet position in z-direction.
    camera_target : :obj:`str` or array_like :obj:`float`, optional
                    Initial camera target. ``auto`` (default) to look towards the system center.
    camera_right : array_like :obj:`float`, optional
                   Camera right vector in system coordinates. Default is [1, 0, 0]
    particle_sizes : :obj:`str` or array_like :obj:`float` or callable, optional
                     auto (default): The Lennard-Jones sigma value of the
                     self-interaction is used for the particle diameter.
                     callable: A lambda function with one argument. Internally,
                     the numerical particle type is passed to the lambda
                     function to determine the particle radius.  list: A list
                     of particle radii, indexed by the particle type.
    particle_coloring : :obj:`str`, optional
                        auto (default): Colors of charged particles are
                        specified by particle_charge_colors, neutral particles
                        by particle_type_colors. charge: Minimum and maximum
                        charge of all particles is determined by the
                        visualizer. All particles are colored by a linear
                        interpolation of the two colors given by
                        particle_charge_colors according to their charge. type:
                        Particle colors are specified by particle_type_colors,
                        indexed by their numerical particle type.
    particle_type_colors : array_like :obj:`float`, optional
                           Colors for particle types.
    particle_type_materials : :obj:`str`, optional
                              Materials of the particle types.
    particle_charge_colors : array_like :obj:`float`, optional
                             Two colors for min/max charged particles.
    draw_constraints : :obj:`bool`, optional
                       Enables constraint visualization. For simple constraints
                       (planes, spheres and cylinders), OpenGL primitives are
                       used. Otherwise, visualization by rasterization is used.
    rasterize_pointsize : :obj:`float`, optional
                          Point size for the rasterization dots.
    rasterize_resolution : :obj:`float`, optional
                           Accuracy of the rasterization.
    quality_constraints : :obj:`int`, optional
                          The number of subdivisions for primitive constraints.
    constraint_type_colors : array_like :obj:`float`, optional
                             Colors of the constaints by type.
    constraint_type_materials : array_like :obj:`str`, optional 
                                Materials of the constraints by type.
    draw_bonds : :obj:`bool`, optional
                 Enables bond visualization.
    bond_type_radius : array_like :obj:`float`, optional
                       Radii of bonds by type.
    bond_type_colors : array_like :obj:`float`, optional
                       Color of bonds by type.
    bond_type_materials : array_like :obj:`float`, optional
                          Materials of bonds by type.
    ext_force_arrows : :obj:`bool`, optional
                       Enables external force visualization.
    ext_force_arrows_scale : :obj:`float`, optional
                             Scale factor of external force arrows for different particle types.
    ext_force_arrows_type_colors : array_like :obj:`float`, optional
                                   Colors of ext_force arrows for different particle types.
    force_arrows : :obj:`bool`, optional
                   Enables particle force visualization.
    force_arrows_scale : :obj:`float`, optional
                         Scale factor of particle force arrows for different particle types.
    force_arrows_type_colors : array_like :obj:`float`, optional
                               Colors of particle force arrows for different particle types.
    velocity_arrows : :obj:`bool`, optional
                       Enables particle velocity visualization.
    velocity_arrows_scale : :obj:`float`, optional
                             Scale factor of particle velocity arrows for different particle types.
    velocity_arrows_type_colors : array_like :obj:`float`, optional
                                  Colors of particle velocity arrows for different particle types.
    drag_enabled : :obj:`bool`, optional
                   Enables mouse-controlled particles dragging (Default: False)
    drag_force : :obj:`bool`, optional
                 Factor for particle dragging
    light_pos : array_like :obj:`float`, optional
                If auto (default) is used, the light is placed dynamically in
                the particle barycenter of the system. Otherwise, a fixed
                coordinate can be set.
    light_colors : array_like :obj:`float`, optional
                   Three lists to specify ambient, diffuse and specular light colors.
    light_brightness : :obj:`float`, optional
                       Brightness (inverse constant attenuation) of the light.
    light_size : :obj:`float`, optional
                 Size (inverse linear attenuation) of the light. If auto
                 (default) is used, the light size will be set to a reasonable
                 value according to the box size at start.
    spotlight_enabled : :obj:`bool`, optional
                        If set to True (default), it enables a spotlight on the
                        camera position pointing in look direction.
    spotlight_colors : array_like :obj:`float`, optional
                       Three lists to specify ambient, diffuse and specular spotlight colors.
    spotlight_angle : :obj:`float`, optional
                      The spread angle of the spotlight in degrees (from 0 to 90).
    spotlight_brightness : :obj:`float`, optional
                           Brightness (inverse constant attenuation) of the spotlight.
    spotlight_focus : :obj:`float`, optional
                      Focus (spot exponent) for the spotlight from 0 (uniform) to 128.

    """

    def __init__(self, system, **kwargs):
        # MATERIALS
        self.materials = {
            'bright':       [[0.9, 0.9, 0.9], [1.0, 1.0, 1.0], [0.8, 0.8, 0.8], 0.6],
            'medium':       [[0.6, 0.6, 0.6], [0.8, 0.8, 0.8], [0.2, 0.2, 0.2], 0.5],
            'dark':         [[0.4, 0.4, 0.4], [0.5, 0.5, 0.5], [0.1, 0.1, 0.1], 0.4],
            'bluerubber':   [[0, 0, 0.05], [0.4, 0.4, 0.5], [0.04, 0.04, 0.7], 0.078125],
            'redrubber':    [[0.05, 0, 0], [0.5, 0.4, 0.4], [0.7, 0.04, 0.04], 0.078125],
            'yellowrubber': [[0.05, 0.05, 0], [0.5, 0.5, 0.4], [0.7, 0.7, 0.04], 0.078125],
            'greenrubber':  [[0, 0.05, 0], [0.4, 0.5, 0.4], [0.04, 0.7, 0.04], 0.078125],
            'whiterubber':  [[0.05, 0.05, 0.05], [0.5, 0.5, 0.5], [0.7, 0.7, 0.7], 0.078125],
            'cyanrubber':   [[0, 0.05, 0.05], [0.4, 0.5, 0.5], [0.04, 0.7, 0.7], 0.078125],
            'blackrubber':  [[0.02, 0.02, 0.02], [0.01, 0.01, 0.01], [0.4, 0.4, 0.4], 0.078125],
            'emerald':      [[0.0215, 0.1745, 0.0215], [0.07568, 0.61424, 0.07568], [0.633, 0.727811, 0.633], 0.6],
            'jade':         [[0.135, 0.2225, 0.1575], [0.54, 0.89, 0.63], [0.316228, 0.316228, 0.316228], 0.1],
            'obsidian':     [[0.05375, 0.05, 0.06625], [0.18275, 0.17, 0.22525], [0.332741, 0.328634, 0.346435], 0.3],
            'pearl':        [[0.25, 0.20725, 0.20725], [1, 0.829, 0.829], [0.296648, 0.296648, 0.296648], 0.088],
            'ruby':         [[0.1745, 0.01175, 0.01175], [0.61424, 0.04136, 0.04136], [0.727811, 0.626959, 0.626959], 0.6],
            'turquoise':    [[0.1, 0.18725, 0.1745], [0.396, 0.74151, 0.69102], [0.297254, 0.30829, 0.306678], 0.1],
            'brass':        [[0.329412, 0.223529, 0.027451], [0.780392, 0.568627, 0.113725], [0.992157, 0.941176, 0.807843], 0.21794872],
            'bronze':       [[0.2125, 0.1275, 0.054], [0.714, 0.4284, 0.18144], [0.393548, 0.271906, 0.166721], 0.2],
            'chrome':       [[0.25, 0.25, 0.25], [0.4, 0.4, 0.4], [0.774597, 0.774597, 0.774597], 0.6],
            'copper':       [[0.19125, 0.0735, 0.0225], [0.7038, 0.27048, 0.0828], [0.256777, 0.137622, 0.086014], 0.1],
            'gold':         [[0.24725, 0.1995, 0.0745], [0.75164, 0.60648, 0.22648], [0.628281, 0.555802, 0.366065], 0.4],
            'silver':       [[0.19225, 0.19225, 0.19225], [0.50754, 0.50754, 0.50754], [0.508273, 0.508273, 0.508273], 0.4],
            'blackplastic': [[0, 0, 0], [0.01, 0.01, 0.01], [0.5, 0.5, 0.5], 0.25],
            'cyanplastic':  [[0, 0.1, 0.06], [0, 0.50980392, 0.50980392], [0.50196078, 0.50196078, 0.50196078], 0.25],
            'greenplastic': [[0, 0, 0], [0.1, 0.35, 0.1], [0.45, 0.55, 0.45], 0.25],
            'redplastic':   [[0, 0, 0], [0.5, 0, 0], [0.7, 0.6, 0.6], 0.25],
            'whiteplastic': [[0, 0, 0], [0.55, 0.55, 0.55], [0.7, 0.7, 0.7], 0.25],
            'yellowplastic': [[0, 0, 0], [0.5, 0.5, 0], [0.6, 0.6, 0.5], 0.25]}

        # DEFAULT PROPERTIES
        self.specs = {
            'window_size': [800, 800],
            'name': 'Espresso Visualization',
            'background_color': [0, 0, 0],
            'update_fps': 30.0,
            'periodic_images': [0, 0, 0],
            'draw_box': True,
            'draw_axis': True,
            'quality_particles': 20,
            'quality_bonds': 16,
            'quality_arrows': 16,
            'arrows_radius': 0.25,
            'quality_constraints': 32,
            'close_cut_distance': 0.1,
            'far_cut_distance': 5,
            'camera_position': 'auto',
            'camera_target': 'auto',
            'camera_right': [1.0, 0.0, 0.0],

            'particle_coloring': 'auto',
            'particle_sizes': 'auto',
            'particle_type_colors': [[1, 1, 0, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            'particle_type_materials': ['medium'],
            'particle_charge_colors': [[1, 0, 0, 1], [0, 1, 0, 1]],

            'draw_constraints': True,
            'rasterize_pointsize': 10,
            'rasterize_resolution': 75.0,
            'constraint_type_colors': [[0.5, 0.5, 0.5, 0.9], [0, 0.5, 0.5, 0.9], [0.5, 0, 0.5, 0.9], [0.5, 0.5, 0, 0.9], [0, 0, 0.5, 0.9], [0.5, 0, 0, 0.9]],
            'constraint_type_materials': ['medium'],

            'draw_bonds': True,
            'bond_type_radius': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            'bond_type_colors': [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            'bond_type_materials': ['medium'],

            'ext_force_arrows': False,
            'ext_force_arrows_scale': [1.0],
            'ext_force_arrows_type_colors': [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            
            'velocity_arrows': False,
            'velocity_arrows_scale': [1.0],
            'velocity_arrows_type_colors': [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            
            'force_arrows': False,
            'force_arrows_scale': [1.0],
            'force_arrows_type_colors': [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],

            'LB': False,
            'LB_plane_axis': 2,
            'LB_plane_dist': 0,
            'LB_plane_ngrid': 5,
            'LB_vel_scale': 1.0,

            'light_pos': 'auto',
            'light_colors': [[0.1, 0.1, 0.2, 1.0], [0.9, 0.9, 0.9, 1.0], [1.0, 1.0, 1.0, 1.0]],
            'light_brightness': 1.0,
            'light_size': 'auto',

            'spotlight_enabled': True,
            'spotlight_colors': [[0.2, 0.2, 0.3, 1.0], [0.5, 0.5, 0.5, 1.0], [1.0, 1.0, 1.0, 1.0]],
            'spotlight_angle': 45,
            'spotlight_focus': 1,
            'spotlight_brightness': 0.6,

            'drag_enabled': False,
            'drag_force': 3.0
        }

        # OVERWRITE WITH USER PROPERTIES
        for key in kwargs:
            if key not in self.specs:
                raise ValueError(
                    key + ' is not a valid visualization property')
            else:
                self.specs[key] = kwargs[key]

        # DEPENDENCIES
        IF not EXTERNAL_FORCES:
            self.specs['drag_enabled'] = False

        self.invBackgroundCol = np.array([1 - self.specs['background_color'][0], 1 -
                                          self.specs['background_color'][1], 1 - self.specs['background_color'][2]])

        self.system = system
        self.started = False
        self.quit_savely = False
        self.keyboardManager = KeyboardManager()
        self.mouseManager = MouseManager()
        self.timers = []

    def register_callback(self, cb, interval=1000):
        """Register timed callbacks.
        """

        self.timers.append((int(interval), cb))

    def run(self, integ_steps=1):
        """Convenience method wiwith a simple integration thread.
        """

        def main():
            while True:
                try:
                    self.system.integrator.run(integ_steps)
                    self.update()
                except:
                    os._exit(1)

        t = Thread(target=main)
        t.daemon = True
        t.start()

        self.start()

    def start(self):
        """The blocking start method.
        """

        self.init_opengl()
        self.init_espresso_visualization()
        self.init_camera()
        self.init_controls()
        self.init_callbacks()

        # POST DISPLAY WITH 60FPS
        def timed_update_redraw(data):
            glutPostRedisplay()
            self.keyboardManager.handle_input()
            glutTimerFunc(17, timed_update_redraw, -1)

        # PLACE LIGHT AT PARTICLE CENTER, DAMPED SPRING FOR SMOOTH POSITION
        # CHANGE, CALL WITH 10FPS
        def timed_update_centerLight(data):
            if self.hasParticleData:
                ldt = 0.8
                cA = (self.particle_COM - self.smooth_light_pos) * \
                    0.1 - self.smooth_light_posV * 1.8
                self.smooth_light_posV += ldt * cA
                self.smooth_light_pos += ldt * self.smooth_light_posV
                self.updateLightPos = True
            glutTimerFunc(100, timed_update_centerLight, -1)

        # AVERAGE PARTICLE COM ONLY EVERY 2sec
        def timed_update_particleCOM(data):
            if self.hasParticleData:
                if len(self.particles['coords']) > 0:
                    self.particle_COM = np.average(
                        self.particles['coords'], axis=0)
            glutTimerFunc(2000, timed_update_particleCOM, -1)

        self.started = True
        self.hasParticleData = False

        glutTimerFunc(17, timed_update_redraw, -1)

        if self.specs['light_pos'] == 'auto':
            glutTimerFunc(2000, timed_update_particleCOM, -1)
            glutTimerFunc(60, timed_update_centerLight, -1)
        # START THE BLOCKING MAIN LOOP
        glutMainLoop()

    def update(self):
        """Update method to be called after integration.
        Changes of espresso system can only happen here.

        """
        if self.started:

            # UPDATE ON STARTUP
            if not self.hasParticleData:
                self.update_particles()
                self.update_charge_color_range()
                self.update_bonds()
                IF CONSTRAINTS:
                    self.update_constraints()
                self.hasParticleData = True

            # IF CALLED TOO OFTEN, ONLY UPDATE WITH GIVEN FREQ
            self.elapsedTime += (time.time() - self.measureTimeBeforeIntegrate)
            if self.elapsedTime > 1.0 / self.specs['update_fps']:
                self.elapsedTime = 0
                self.update_particles()
                if self.specs['LB']:
                    self.update_lb()
                # KEYBOARD CALLBACKS MAY CHANGE ESPRESSO SYSTEM PROPERTIES,
                # ONLY SAVE TO CHANGE HERE
                for c in self.keyboardManager.userCallbackStack:
                    c()
                self.keyboardManager.userCallbackStack = []

            self.measureTimeBeforeIntegrate = time.time()

            IF EXTERNAL_FORCES:
                if self.triggerSetParticleDrag == True and self.dragId != -1:
                    self.system.part[self.dragId].ext_force = self.dragExtForce
                    self.triggerSetParticleDrag = False
                elif self.triggerResetParticleDrag == True and self.dragId != -1:
                    self.system.part[self.dragId].ext_force = self.extForceOld
                    self.triggerResetParticleDrag = False
                    self.dragId = -1

        # Escape was pressed: wait for ES to finish, then call sys exit from
        # main thread
        if self.quit_savely:
            os._exit(1)

    # GET THE PARTICLE DATA
    def update_particles(self):
        IF EXTERNAL_FORCES and ELECTROSTATICS:
            self.particles = {'coords': self.system.part[:].pos_folded,
                              'velocities': self.system.part[:].v if self.specs['velocity_arrows'] else [0, 0, 0] * len(self.system.part),
                              'forces': self.system.part[:].f if self.specs['force_arrows'] else [0, 0, 0] * len(self.system.part),
                              'types': self.system.part[:].type,
                              'ext_forces': self.system.part[:].ext_force,
                              'charges': self.system.part[:].q}
        ELIF EXTERNAL_FORCES and not ELECTROSTATICS:
            self.particles = {'coords': self.system.part[:].pos_folded,
                              'velocities': self.system.part[:].v if self.specs['velocity_arrows'] else [0, 0, 0] * len(self.system.part),
                              'forces': self.system.part[:].f if self.specs['force_arrows'] else [0, 0, 0] * len(self.system.part),
                              'types': self.system.part[:].type,
                              'ext_forces': self.system.part[:].ext_force,
                              'charges': [0] * len(self.system.part)}
        ELIF not EXTERNAL_FORCES and ELECTROSTATICS:
            self.particles = {'coords': self.system.part[:].pos_folded,
                              'velocities': self.system.part[:].v if self.specs['velocity_arrows'] else [0, 0, 0] * len(self.system.part),
                              'forces': self.system.part[:].f if self.specs['force_arrows'] else [0, 0, 0] * len(self.system.part),
                              'types': self.system.part[:].type,
                              'ext_forces': [0, 0, 0] * len(self.system.part),
                              'charges': self.system.part[:].q}
        ELIF not EXTERNAL_FORCES and not ELECTROSTATICS:
            self.particles = {'coords': self.system.part[:].pos_folded,
                              'velocities': self.system.part[:].v if self.specs['velocity_arrows'] else [0, 0, 0] * len(self.system.part),
                              'forces': self.system.part[:].f if self.specs['force_arrows'] else [0, 0, 0] * len(self.system.part),
                              'types': self.system.part[:].type,
                              'ext_forces': [0, 0, 0] * len(self.system.part),
                              'charges': [0] * len(self.system.part)}

    def update_lb(self):
        agrid = self.lb_params['agrid']
        self.lb_plane_vel = []
        ng = self.specs['LB_plane_ngrid']
        for xi in xrange(ng):
            for xj in xrange(ng):
                pp = np.copy((self.lb_plane_p + xi * 1.0 / ng * self.lb_plane_b1 +
                      xj * 1.0 / ng * self.lb_plane_b2) % self.system.box_l)
                i, j, k = (int(ppp / agrid) for ppp in pp)
                lb_vel = np.copy(self.lb[i, j, k].velocity)
                self.lb_plane_vel.append([pp, lb_vel])

    def edges_from_pn(self, p, n, diag):
        v1, v2 = self.get_tangents(n)

        edges = []
        edges.append(p + diag * v1)
        edges.append(p + diag * v2)
        edges.append(p - diag * v1)
        edges.append(p - diag * v2)
        return edges

    # GET THE update_constraints DATA
    def update_constraints(self):

        box_diag = pow(pow(self.system.box_l[0], 2) + pow(
            self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)

        self.shapes = collections.defaultdict(list)

        # Collect shapes and interaction type (for coloring) from constraints
        coll_shape_obj = collections.defaultdict(list)
        for c in self.system.constraints:
            if type(c) == espressomd.constraints.ShapeBasedConstraint:
                t = c.get_parameter('particle_type')
                s = c.get_parameter('shape')
                n = s.name()
                if n in ['Shapes::Wall', 'Shapes::Cylinder', 'Shapes::Sphere', 'Shapes::SpheroCylinder']:
                    coll_shape_obj[n].append([s, t])
                else:
                    coll_shape_obj['Shapes::Misc'].append([s, t])

        # TODO: get shapes from lbboundaries
        for s in coll_shape_obj['Shapes::Wall']:
            d = s[0].get_parameter('dist')
            n = s[0].get_parameter('normal')
            edges = self.edges_from_pn(d * np.array(n), n, 2 * box_diag)
            self.shapes['Shapes::Wall'].append([edges, s[1]])

        for s in coll_shape_obj['Shapes::Cylinder']:
            pos = np.array(s[0].get_parameter('center'))
            a = np.array(s[0].get_parameter('axis'))
            l = s[0].get_parameter('length')
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::Cylinder'].append(
                [pos - a * l * 0.5, pos + a * l * 0.5, r, s[1]])

        for s in coll_shape_obj['Shapes::Sphere']:
            pos = np.array(s[0].get_parameter('center'))
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::Sphere'].append([pos, r, s[1]])

        for s in coll_shape_obj['Shapes::SpheroCylinder']:
            pos = np.array(s[0].get_parameter('center'))
            a = np.array(s[0].get_parameter('axis'))
            l = s[0].get_parameter('length')
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::SpheroCylinder'].append(
                [pos - a * l * 0.5, pos + a * l * 0.5, r, s[1]])

        for s in coll_shape_obj['Shapes::Misc']:
            self.shapes['Shapes::Misc'].append(
                [self.rasterize_brute_force(s[0]), s[1]])

    def get_tangents(self, n):
        n = np.array(n)
        v1 = np.random.randn(3)
        v1 -= v1.dot(n) * n / np.linalg.norm(n)**2
        v2 = np.cross(n, v1)
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        return v1, v2

    def rasterize_brute_force(self, shape):
        sp = max(self.system.box_l) / self.specs['rasterize_resolution']
        res = np.array(self.system.box_l) / sp

        points = []
        for i in range(int(res[0])):
            for j in range(int(res[1])):
                for k in range(int(res[2])):
                    p = np.array([i, j, k]) * sp
                    dist, vec = shape.call_method(
                        "calc_distance", position=p.tolist())
                    if not np.isnan(vec).any() and not np.isnan(dist) and abs(dist) < sp:
                        points.append((p - vec).tolist())
        return points

    # GET THE BOND DATA, SO FAR CALLED ONCE UPON INITIALIZATION
    def update_bonds(self):
        if self.specs['draw_bonds']:
            self.bonds = []
            for i in range(len(self.system.part)):

                bs = self.system.part[i].bonds
                for b in bs:
                    t = b[0].type_number()
                    # b[0]: Bond, b[1:] Partners
                    for p in b[1:]:
                        self.bonds.append([i, p, t])

    # DRAW CALLED AUTOMATICALLY FROM GLUT DISPLAY FUNC
    def draw(self):

        if self.specs['LB']:
            self.draw_lb_vel()
        if self.specs['draw_box']:
            self.draw_system_box()

        if self.specs['draw_axis']:
            axis_fac = 0.2
            axis_r = np.min(self.system.box_l) / 50.0
            draw_arrow([0, 0, 0], [self.system.box_l[0] * axis_fac, 0, 0], axis_r, [
                       1, 0, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, self.system.box_l[2] * axis_fac, 0], axis_r, [
                       0, 1, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, 0, self.system.box_l[2] * axis_fac], axis_r, [
                       0, 0, 1, 1], self.materials['chrome'], self.specs['quality_arrows'])

        self.draw_system_particles()

        if self.specs['draw_bonds']:
            self.draw_bonds()

        IF CONSTRAINTS:
            if self.specs['draw_constraints']:
                self.draw_constraints()

    def draw_system_box(self):
        draw_box([0, 0, 0], self.system.box_l, self.invBackgroundCol)

    def draw_constraints(self):

        for i in range(6):
            glEnable(GL_CLIP_PLANE0 + i)
            glClipPlane(GL_CLIP_PLANE0 + i, self.box_eqn[i])

        for s in self.shapes['Shapes::Sphere']:
            draw_sphere(s[0], s[1], self.modulo_indexing(self.specs['constraint_type_colors'], s[2]), self.materials[self.modulo_indexing(
                self.specs['constraint_type_materials'], s[2])], self.specs['quality_constraints'])

        for s in self.shapes['Shapes::SpheroCylinder']:
            draw_sphero_cylinder(
                s[0], s[1], s[2], self.modulo_indexing(
                    self.specs['constraint_type_colors'], s[3]),
                               self.materials[self.modulo_indexing(self.specs['constraint_type_materials'], s[3])], self.specs['quality_constraints'])

        for s in self.shapes['Shapes::Wall']:
            draw_plane(
                s[0], self.modulo_indexing(
                    self.specs['constraint_type_colors'], s[1]),
                      self.materials[self.modulo_indexing(self.specs['constraint_type_materials'], s[1])])

        for s in self.shapes['Shapes::Cylinder']:
            draw_cylinder(s[0], s[1], s[2], self.modulo_indexing(self.specs['constraint_type_colors'], s[3]), self.materials[self.modulo_indexing(
                self.specs['constraint_type_materials'], s[3])], self.specs['quality_constraints'], True)

        for s in self.shapes['Shapes::Misc']:
            draw_points(s[0], self.specs['rasterize_pointsize'],  self.modulo_indexing(
                self.specs['constraint_type_colors'], s[1]), self.materials[self.modulo_indexing(self.specs['constraint_type_materials'], s[1])])

        for i in range(6):
            glDisable(GL_CLIP_PLANE0 + i)

    def determine_radius(self, ptype):
        def radiusByLJ(ptype):
            try:
                radius = self.system.non_bonded_inter[ptype, ptype].lennard_jones.get_params()[
                    'sigma'] * 0.5
            except BaseException:
                radius = 0.5
            if radius == 0:
                radius = 0.5
            return radius

        if self.specs['particle_sizes'] == 'auto':
            radius = radiusByLJ(ptype)
        elif callable(self.specs['particle_sizes']):
            radius = self.specs['particle_sizes'](ptype)
        else:
            try:
                radius = self.modulo_indexing(
                    self.specs['particle_sizes'], ptype)
            except BaseException:
                radius = self.radiusByLJ(ptype)
        return radius

    def draw_system_particles(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        for pid in pIds:
            pos = coords[pid]
            q = self.particles['charges'][pid]
            ptype = int(self.particles['types'][pid])
            ext_f = self.particles['ext_forces'][pid]

            radius = self.determine_radius(ptype)

            m = self.modulo_indexing(
                self.specs['particle_type_materials'], ptype)
            material = self.materials[m]

            if self.specs['particle_coloring'] == 'id':
                color = self.id_to_fcolor(pid)
                glColor(color)
            elif self.specs['particle_coloring'] == 'auto':
                # Color auto: Charge then Type
                if q != 0:
                    color = self.color_by_charge(q)
                else:
                    color = self.modulo_indexing(
                        self.specs['particle_type_colors'], ptype)
            elif self.specs['particle_coloring'] == 'charge':
                color = self.color_by_charge(q)
            elif self.specs['particle_coloring'] == 'type':
                color = self.modulo_indexing(
                    self.specs['particle_type_colors'], ptype)

            draw_sphere(pos, radius, color, material,
                        self.specs['quality_particles'])
            for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                    for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                        if imx != 0 or imy != 0 or imz != 0:
                            redraw_sphere(
                                pos + (imx * self.imPos[0] + imy * self.imPos[1] + imz * self.imPos[2]), radius, self.specs['quality_particles'])

            IF EXTERNAL_FORCES:
                if self.specs['ext_force_arrows'] or pid == self.dragId:
                    if ext_f[0] != 0 or ext_f[1] != 0 or ext_f[2] != 0:
                        if pid == self.dragId:
                            sc = 1
                        else:
                            sc = self.modulo_indexing(
                                self.specs['ext_force_arrows_scale'], ptype)
                        if sc > 0:
                            col = self.modulo_indexing(self.specs['ext_force_arrows_type_colors'], ptype)
                            draw_arrow(pos, np.array(ext_f) * sc, self.specs['arrows_radius']
                                       , col, self.materials['chrome'], self.specs['quality_arrows'])
            
            if self.specs['velocity_arrows']:
                sc = self.modulo_indexing(self.specs['velocity_arrows_scale'], ptype)
                if sc > 0:
                    v = self.particles['velocities'][pid]
                    col = self.modulo_indexing(self.specs['velocity_arrows_type_colors'], ptype)
                    draw_arrow(pos, np.array(v) * sc, self.specs['arrows_radius']
                               , col, self.materials['chrome'], self.specs['quality_arrows'])
            
            if self.specs['force_arrows']:
                sc = self.modulo_indexing(self.specs['force_arrows_scale'], ptype)
                if sc > 0:
                    v = self.particles['forces'][pid]
                    col = self.modulo_indexing(self.specs['force_arrows_type_colors'], ptype)
                    draw_arrow(pos, np.array(v) * sc, self.specs['arrows_radius']
                               , col, self.materials['chrome'], self.specs['quality_arrows'])
                

    def draw_bonds(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        b2 = self.system.box_l[0] / 2.0
        box_l2_sqr = pow(b2, 2.0)
        for b in self.bonds:
            col = self.modulo_indexing(self.specs['bond_type_colors'], b[2])
            mat = self.materials[self.modulo_indexing(
                self.specs['bond_type_materials'], b[2])]
            radius = self.modulo_indexing(
                self.specs['bond_type_radius'], b[2])
            d = coords[b[0]] - coords[b[1]]
            bondLen_sqr = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]

            if bondLen_sqr < box_l2_sqr:
                draw_cylinder(coords[b[0]], coords[b[1]], radius,
                              col, mat, self.specs['quality_bonds'])
                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                im = np.array([imx, imy, imz])
                                draw_cylinder(coords[b[0]] + im * self.imPos[dim], coords[b[1]] +
                                              im * self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])
            else:
                l = coords[b[0]] - coords[b[1]]
                l0 = coords[b[0]]
                hits = 0
                for i in range(6):
                    lineBoxNDot = float(np.dot(l, self.box_n[i]))
                    if lineBoxNDot == 0:
                        continue
                    s = l0 - \
                        np.dot(l0 - self.box_p[i],
                               self.box_n[i]) / lineBoxNDot * l
                    if self.is_inside_box(s):
                        if lineBoxNDot < 0:
                            s0 = s
                        else:
                            s1 = s
                        hits += 1
                        if hits >= 2:
                            break
                draw_cylinder(coords[b[0]], s0, radius, col,
                              mat, self.specs['quality_bonds'])
                draw_cylinder(coords[b[1]], s1, radius, col,
                              mat, self.specs['quality_bonds'])

                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                im = np.array([imx, imy, imz])
                                draw_cylinder(coords[b[0]] + im * self.imPos[dim], s0 + im *
                                              self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])
                                draw_cylinder(coords[b[1]] + im * self.imPos[dim], s1 + im *
                                              self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])

    # HELPER TO DRAW PERIODIC BONDS
    def is_inside_box(self, p):
        eps = 1e-5
        for i in range(3):
            if p[i] < -eps or p[i] > eps + self.system.box_l[i]:
                return False
        return True

# VOXELS FOR LB VELOCITIES
    def draw_lb_vel(self):

        for lbl in self.lb_plane_vel:
            p = lbl[0]
            v = lbl[1]
            c = np.linalg.norm(v)
            draw_arrow(p, v * self.specs['LB_vel_scale'], self.lb_arrow_radius, [
                       1, 1, 1, 1], self.materials['chrome'], 16)

    # USE MODULO IF THERE ARE MORE PARTICLE TYPES THAN TYPE DEFINITIONS FOR
    # COLORS, MATERIALS ETC..
    def modulo_indexing(self, l, t):
        return l[t % len(l)]

    # FADE PARTICE CHARGE COLOR FROM WHITE (q=0) to PLUSCOLOR (q=q_max) RESP
    # MINUSCOLOR (q=q_min)
    def color_by_charge(self, q):
        if q < 0:
            c = 1.0 * q / self.minq
            return np.array(self.specs['particle_charge_colors'][0]) * c + (1 - c) * np.array([1, 1, 1, 1])
        else:
            c = 1.0 * q / self.maxq

            return np.array(self.specs['particle_charge_colors'][1]) * c + (1 - c) * np.array([1, 1, 1, 1])

    # ON INITIALIZATION, CHECK q_max/q_min
    def update_charge_color_range(self):
        if len(self.particles['charges'][:]) > 0:
            self.minq = min(self.particles['charges'][:])
            self.maxq = max(self.particles['charges'][:])

    # INITS FOR GLUT FUNCTIONS
    def init_callbacks(self):
        # OpenGl Callbacks
        def display():
            if self.hasParticleData:
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

                glLoadMatrixf(self.camera.modelview)

                if self.updateLightPos:
                    self.set_light_pos()
                    self.updateLightPos = False

                self.draw()

                glutSwapBuffers()
            return

        def keyboard_up(button, x, y):
            if type(button) is bytes:
                button = button.decode("utf-8")
            self.keyboardManager.keyboard_up(button)
            return

        def keyboard_down(button, x, y):
            if type(button) is bytes:
                button = button.decode("utf-8")
            self.keyboardManager.keyboard_down(button)
            return

        def mouse(button, state, x, y):
            self.mouseManager.mouse_click(button, state, x, y)
            return

        def motion(x, y):
            self.mouseManager.mouse_move(x, y)
            return

        def idle_update():
            return

        # CALLED ION WINDOW POSITION/SIZE CHANGE
        def reshape_window(w, h):
            glViewport(0, 0, w, h)
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            box_diag = pow(pow(self.system.box_l[0], 2) + pow(
                self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)
            gluPerspective(
                40, 1.0 * w / h, self.specs['close_cut_distance'], self.specs['far_cut_distance'] * box_diag)
            glMatrixMode(GL_MODELVIEW)
            self.specs['window_size'][0] = 1.0 * w
            self.specs['window_size'][1] = 1.0 * h
            # glPushMatrix()

        def close_window():
            os._exit(1)

        # TIMERS FOR register_callback
        def dummy_timer(index):
            self.timers[index][1]()
            glutTimerFunc(self.timers[index][0], dummy_timer, index)

        glutDisplayFunc(display)
        glutMouseFunc(mouse)
        glutKeyboardFunc(keyboard_down)
        glutKeyboardUpFunc(keyboard_up)
        glutReshapeFunc(reshape_window)
        glutMotionFunc(motion)
        glutWMCloseFunc(close_window)

        glutIdleFunc(idle_update)

        index = 0
        for t in self.timers:
            glutTimerFunc(t[0], dummy_timer, index)
            index += 1

    # CLICKED ON PARTICLE: DRAG; CLICKED ON BACKGROUND: CAMERA
    def mouse_motion(self, mousePos, mousePosOld, mouseButtonState):

        if self.dragId != -1:
            ppos = self.particles['coords'][self.dragId]
            viewport = glGetIntegerv(GL_VIEWPORT)
            mouseWorld = gluUnProject(
                mousePos[0], viewport[3] - mousePos[1], self.depth)

            self.dragExtForce = self.specs['drag_force'] * \
                (np.asarray(mouseWorld) - np.array(ppos))
            self.triggerSetParticleDrag = True
        else:
            self.camera.rotate_camera(mousePos, mousePosOld, mouseButtonState)

    # DRAW SCENE AGAIN WITHOUT LIGHT TO IDENTIFY PARTICLE ID BY PIXEL COLOR
    def set_particle_drag(self, pos, pos_old):

        glClearColor(0.0, 0.0, 0.0, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glLoadMatrixf(self.camera.modelview)

        oldColMode = self.specs['particle_coloring']
        self.specs['particle_coloring'] = 'id'
        glDisable(GL_LIGHTING)
        glDisable(GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            glDisable(GL_LIGHT1)
        self.draw_system_particles()
        viewport = glGetIntegerv(GL_VIEWPORT)

        readPixel = glReadPixelsui(
            pos[0], viewport[3] - pos[1], 1, 1, GL_RGB, GL_FLOAT)[0][0]
        depth = glReadPixelsf(
            pos[0], viewport[3] - pos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]
        pid = self.fcolor_to_id(readPixel)

        self.dragId = pid
        if pid != -1:
            self.dragPosInitial = self.particles['coords'][self.dragId]
            self.extForceOld = self.particles['ext_forces'][self.dragId][:]
            self.depth = depth
        self.specs['particle_coloring'] = oldColMode
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            glEnable(GL_LIGHT1)
        glClearColor(self.specs['background_color'][0],
                     self.specs['background_color'][1],
                     self.specs['background_color'][2], 1.)

    def reset_particle_drag(self, pos, pos_old):
        if self.dragId != -1:
            self.triggerResetParticleDrag = True

    def id_to_fcolor(self, pid):
        pid += 1
        return [int(pid / (256 * 256)) / 255.0, int((pid % (256 * 256)) / 256) / 255.0, (pid % 256) / 255.0, 1.0]

    def fcolor_to_id(self, fcol):
        if (fcol == [0, 0, 0]).all():
            return -1
        else:
            return 256 * 256 * int(fcol[0] * 255) + 256 * int(fcol[1] * 255) + int(fcol[2] * 255) - 1

    # ALL THE INITS
    def init_espresso_visualization(self):
        self.maxq = 0
        self.minq = 0

        self.dragId = -1
        self.dragPosInitial = []
        self.extForceOld = []
        self.dragExtForceOld = []
        self.triggerResetParticleDrag = False
        self.triggerSetParticleDrag = False

        self.depth = 0

        self.imPos = [np.array([self.system.box_l[0], 0, 0]), np.array(
            [0, self.system.box_l[1], 0]), np.array([0, 0, self.system.box_l[2]])]

        if self.specs['LB']:
            for a in self.system.actors:
                pa = a.get_params()
                if 'agrid' in pa:
                    self.lb_params = pa
                    self.lb = a
                    break

            if self.specs['LB_plane_axis'] == 0:
                pn = [1.0, 0.0, 0.0]
                self.lb_plane_b1 = [0.0,1.0,0.0]
                self.lb_plane_b2 = [0.0,0.0,1.0]
            elif self.specs['LB_plane_axis'] == 1:
                pn = [0.0, 1.0, 0.0]
                self.lb_plane_b1 = [1.0,0.0,0.0]
                self.lb_plane_b2 = [0.0,0.0,1.0]
            else:
                pn = [0.0, 0.0, 1.0]
                self.lb_plane_b1 = [1.0,0.0,0.0]
                self.lb_plane_b2 = [0.0,1.0,0.0]

            self.lb_plane_b1 *= self.system.box_l
            self.lb_plane_b2 *= self.system.box_l
            self.lb_plane_p = np.array(pn) * self.specs['LB_plane_dist']

            self.lb_arrow_radius = self.system.box_l[
                self.specs['LB_plane_axis']] * 0.005

            self.lb_min_vel = np.array([-1e-6] * 3)
            self.lb_max_vel = np.array([1e-6] * 3)
            self.lb_vel_range = self.lb_max_vel - self.lb_min_vel
            self.lb_min_dens = np.array([0] * 3)
            self.lb_max_dens = np.array([0] * 3)

            self.update_lb()

        self.elapsedTime = 0
        self.measureTimeBeforeIntegrate = 0

        self.box_size_dependence()

    # BOX PLANES (NORMAL, ORIGIN) FOR PERIODIC BONDS
    def box_size_dependence(self):
        self.box_n = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array(
            [0, 0, 1]), np.array([-1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, -1])]
        self.box_p = [np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array(
            self.system.box_l), np.array(self.system.box_l), np.array(self.system.box_l)]
        self.box_eqn = []
        self.box_eqn.append(
            (self.box_n[0][0], self.box_n[0][1], self.box_n[0][2], self.system.box_l[0] * 0.001))
        self.box_eqn.append(
            (self.box_n[1][0], self.box_n[1][1], self.box_n[1][2], self.system.box_l[1] * 0.001))
        self.box_eqn.append(
            (self.box_n[2][0], self.box_n[2][1], self.box_n[2][2], self.system.box_l[2] * 0.001))
        self.box_eqn.append(
            (self.box_n[3][0], self.box_n[3][1], self.box_n[3][2], self.system.box_l[0] * 1.001))
        self.box_eqn.append(
            (self.box_n[4][0], self.box_n[4][1], self.box_n[4][2], self.system.box_l[1] * 1.001))
        self.box_eqn.append(
            (self.box_n[5][0], self.box_n[5][1], self.box_n[5][2], self.system.box_l[2] * 1.001))

    # DEFAULT CONTROLS
    def init_controls(self):
        # MOUSE LOOK/ROTATE/DRAG
        self.mouseManager.register_button(MouseButtonEvent(
            None, MouseFireEvent.FreeMotion, self.mouse_motion))

        self.mouseManager.register_button(MouseButtonEvent(
            3, MouseFireEvent.ButtonPressed, self.camera.move_backward))

        self.mouseManager.register_button(MouseButtonEvent(
            4, MouseFireEvent.ButtonPressed, self.camera.move_forward))

        # START/STOP DRAG
        if self.specs['drag_enabled']:
            self.mouseManager.register_button(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonPressed, self.set_particle_drag, True))
            self.mouseManager.register_button(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonReleased, self.reset_particle_drag, True))

        # KEYBOARD BUTTONS
        self.keyboardManager.register_button(KeyboardButtonEvent(
            '\x1b', KeyboardFireEvent.Pressed, self.quit, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'w', KeyboardFireEvent.Hold, self.camera.move_up, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            's', KeyboardFireEvent.Hold, self.camera.move_down, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'a', KeyboardFireEvent.Hold, self.camera.move_left, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'd', KeyboardFireEvent.Hold, self.camera.move_right, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'e', KeyboardFireEvent.Hold, self.camera.rotate_system_XR, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'q', KeyboardFireEvent.Hold, self.camera.rotate_system_XL, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'z', KeyboardFireEvent.Hold, self.camera.rotate_system_YR, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'c', KeyboardFireEvent.Hold, self.camera.rotate_system_YL, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'r', KeyboardFireEvent.Hold, self.camera.rotate_system_ZR, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'f', KeyboardFireEvent.Hold, self.camera.rotate_system_ZL, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            't', KeyboardFireEvent.Hold, self.camera.move_forward, True))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            'g', KeyboardFireEvent.Hold, self.camera.move_backward, True))

    # CALLED ON ESCAPE PRESSED. TRIGGERS sys.exit() after ES is done
    def quit(self):
        self.quit_savely = True

    # ASYNCHRONOUS PARALLEL CALLS OF glLight CAUSES SEG FAULTS, SO ONLY CHANGE
    # LIGHT AT CENTRAL display METHOD AND TRIGGER CHANGES
    def set_light_pos(self):
        if self.specs['light_pos'] == 'auto':
            glLightfv(GL_LIGHT0, GL_POSITION, [
                      self.smooth_light_pos[0], self.smooth_light_pos[1], self.smooth_light_pos[2], 0.6])

        self.set_camera_spotlight()

    def set_camera_spotlight(self):
        # p = np.linalg.norm(self.camera.state_pos) * self.camera.state_target
        p = self.camera.camPos
        fp = [p[0], p[1], p[2], 1]
        glLightfv(GL_LIGHT1, GL_POSITION, fp)
        glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, self.camera.state_target)

    def trigger_light_pos_update(self):
        self.updateLightPos = True

    def init_camera(self):
        b = np.array(self.system.box_l)
        box_diag = np.linalg.norm(b)
        box_center = b * 0.5
        if self.specs['camera_position'] == 'auto':
            cp = [box_center[0], box_center[1], b[2] * 3]
        else:
            cp = self.specs['camera_position']

        if self.specs['camera_target'] == 'auto':
            ct = box_center
        else:
            ct = self.specs['camera_target']

        cr = np.array(self.specs['camera_right'])

        self.camera = Camera(camPos=np.array(cp), camTarget=ct, camRight=cr, moveSpeed=0.5 *
                              box_diag / 17.0,  center=box_center, updateLights=self.trigger_light_pos_update)
        self.smooth_light_pos = np.copy(box_center)
        self.smooth_light_posV = np.array([0.0, 0.0, 0.0])
        self.particle_COM = np.copy(box_center)
        self.set_camera_spotlight()
        self.updateLightPos = True

    def init_opengl(self):
        glutInit(self.specs['name'])
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(self.specs['window_size'][
                           0], self.specs['window_size'][1])

        # glutCreateWindow(bytes(self.specs['name'], encoding='ascii'))
        glutCreateWindow(b"ESPResSo visualization")

        glClearColor(self.specs['background_color'][0], self.specs[
                     'background_color'][1], self.specs['background_color'][2], 1.)

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        glEnable(GL_BLEND)

        # glEnable(GL_CULL_FACE)

        glLineWidth(2.0)
        glutIgnoreKeyRepeat(1)

        # setup lighting
        if self.specs['light_size'] == 'auto':
            box_diag = pow(pow(self.system.box_l[0], 2) + pow(
                self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)
            self.specs['light_size'] = box_diag * 2.0

        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)

        # LIGHT0
        if self.specs['light_pos'] != 'auto':
            glLightfv(
                GL_LIGHT0, GL_POSITION, np.array(
                    self.specs['light_pos']).tolist())
        else:
            glLightfv(
                GL_LIGHT0,
                GL_POSITION,
                (np.array(
                    self.system.box_l) *
                    0.5).tolist())

        glLightfv(GL_LIGHT0, GL_AMBIENT, self.specs['light_colors'][0])
        glLightfv(GL_LIGHT0, GL_DIFFUSE, self.specs['light_colors'][1])
        glLightfv(GL_LIGHT0, GL_SPECULAR, self.specs['light_colors'][2])

        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,
                 1.0 / self.specs['light_brightness'])
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION,
                 1.0 / self.specs['light_size'])
        glEnable(GL_LIGHT0)

        # LIGHT1: SPOTLIGHT ON CAMERA IN LOOK DIRECTION
        if self.specs['spotlight_enabled']:
            glLightfv(GL_LIGHT1, GL_POSITION, [0, 0, 0, 1])

            glLightfv(GL_LIGHT1, GL_AMBIENT, self.specs['spotlight_colors'][0])
            glLightfv(GL_LIGHT1, GL_DIFFUSE, self.specs['spotlight_colors'][1])
            glLightfv(GL_LIGHT1, GL_SPECULAR,
                      self.specs['spotlight_colors'][2])

            glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, self.specs['spotlight_angle'])
            glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, [1.0, 1.0, 1.0])
            glLightf(GL_LIGHT1, GL_SPOT_EXPONENT,
                     self.specs['spotlight_focus'])

            glLightf(GL_LIGHT1, GL_CONSTANT_ATTENUATION,
                     1.0 / self.specs['spotlight_brightness'])
            glLightf(GL_LIGHT1, GL_LINEAR_ATTENUATION, 0.0)
            glLightf(GL_LIGHT1, GL_QUADRATIC_ATTENUATION, 0.0)
            glEnable(GL_LIGHT1)

# END OF MAIN CLASS

# OPENGL DRAW WRAPPERS


def set_solid_material(r, g, b, a=1.0, ambient=[0.6, 0.6, 0.6], diffuse=[1.0, 1.0, 1.0], specular=[0.1, 0.1, 0.1], shininess=0.4):
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  [
                 ambient[0] * r, ambient[1] * g, ambient[2] * g, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [
                 diffuse[0] * r, diffuse[1] * g, diffuse[2] * b, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [
                 specular[0] * r, specular[1] * g, specular[2] * g, a])
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, int(shininess * 128))


def draw_box(p0, s, color):
    set_solid_material(color[0], color[1], color[2])
    glPushMatrix()
    glTranslatef(p0[0], p0[1], p0[2])
    glBegin(GL_LINE_LOOP)
    glVertex3f(0.0, 0.0, 0.0)
    glVertex3f(s[0], 0.0, 0.0)
    glVertex3f(s[0], s[1], 0.0)
    glVertex3f(0, s[1], 0.0)
    glEnd()
    glBegin(GL_LINE_LOOP)
    glVertex3f(0.0, 0.0, s[2])
    glVertex3f(s[0], 0.0, s[2])
    glVertex3f(s[0], s[1], s[2])
    glVertex3f(0, s[1], s[2])
    glEnd()
    glBegin(GL_LINES)
    glVertex3f(0.0, 0.0, 0.0)
    glVertex3f(0.0, 0.0, s[2])
    glVertex3f(s[0], 0.0, 0.0)
    glVertex3f(s[0], 0.0, s[2])
    glVertex3f(s[0], s[1], 0.0)
    glVertex3f(s[0], s[1], s[2])
    glVertex3f(0.0, s[1], 0.0)
    glVertex3f(0.0, s[1], s[2])
    glEnd()
    # glutWireCube(size)
    glPopMatrix()


def draw_sphere(pos, radius, color, material, quality):
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    set_solid_material(color[0], color[1], color[2], color[
                      3], material[0], material[1], material[2], material[3])
    glutSolidSphere(radius, quality, quality)
    glPopMatrix()


def redraw_sphere(pos, radius, quality):
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glutSolidSphere(radius, quality, quality)
    glPopMatrix()


def draw_plane(edges, color, material):

    set_solid_material(color[0], color[1], color[2], color[
        3], material[0], material[1], material[2], material[3])

    glBegin(GL_QUADS)
    for e in edges:
        glVertex3f(e[0], e[1], e[2])
    glEnd()


def draw_cube(pos, size, color, alpha):
    set_solid_material(color[0], color[1], color[2], alpha)
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glutSolidCube(size)
    glPopMatrix()


def draw_triangles(triangles, color, material):
    np.random.seed(1)

    glBegin(GL_TRIANGLES)
    for t in triangles:
        color = np.random.random(3).tolist()
        color.append(1)
        set_solid_material(color[0], color[1], color[2],
                          color[3], material[0], material[1], material[2], material[3])
        for p in t:
            glVertex3f(p[0], p[1], p[2])
    glEnd()


def draw_points(points, pointsize, color, material):
    set_solid_material(color[0], color[1], color[2], color[3],
                      material[0], material[1], material[2], material[3])
    glEnable(GL_POINT_SMOOTH)

    glPointSize(pointsize)
    glBegin(GL_POINTS)
    for p in points:
        glVertex3f(p[0], p[1], p[2])
    glEnd()


def draw_cylinder(posA, posB, radius, color, material, quality, draw_caps=False):
    set_solid_material(color[0], color[1], color[2], color[3],
                      material[0], material[1], material[2], material[3])
    glPushMatrix()
    quadric = gluNewQuadric()

    d = posB - posA

    # angle,t,length = calcAngle(d)
    length = np.linalg.norm(d)
    glTranslatef(posA[0], posA[1], posA[2])

    ax, rx, ry = rotation_helper(d)
    glRotatef(ax, rx, ry, 0.0)
    gluCylinder(quadric, radius, radius, length, quality, quality)

    if draw_caps:
        gluDisk(quadric, 0, radius, quality, quality)
        glTranslatef(0, 0, length)
        gluDisk(quadric, 0, radius, quality, quality)

    glPopMatrix()


def rotation_helper(d):
    if d[2] == 0.0:
        d[2] = 0.0001

    v = np.linalg.norm(d)
    if v == 0:
        ax = 57.2957795
    else:
        ax = 57.2957795 * acos(d[2] / v)

    if d[2] < 0.0:
        ax = -ax
    rx = -d[1] * d[2]
    ry = d[0] * d[2]

    return ax, rx, ry


def draw_sphero_cylinder(posA, posB, radius, color, material, quality):
    set_solid_material(color[0], color[1], color[
        2], color[3], material[0], material[1], material[2], material[3])
    glPushMatrix()
    quadric = gluNewQuadric()

    d = posB - posA
    if d[2] == 0.0:
        d[2] = 0.0001

    v = np.linalg.norm(d)
    if v == 0:
        ax = 57.2957795
    else:
        ax = 57.2957795 * acos(d[2] / v)

    if d[2] < 0.0:
        ax = -ax
    rx = -d[1] * d[2]
    ry = d[0] * d[2]
    length = np.linalg.norm(d)
    glTranslatef(posA[0], posA[1], posA[2])
    glRotatef(ax, rx, ry, 0.0)

    # First hemispherical cap
    clip_plane = GL_CLIP_PLANE0+6
    glEnable(clip_plane)
    glClipPlane(clip_plane, (0, 0, -1, 0))
    gluSphere(quadric, radius, quality, quality)
    glDisable(clip_plane)
    # Cylinder
    gluCylinder(quadric, radius, radius, length, quality, quality)
    # Second hemispherical cap
    glTranslatef(0, 0, v)
    glEnable(clip_plane)
    glClipPlane(clip_plane, (0, 0, 1, 0))
    gluSphere(quadric, radius, quality, quality)
    glDisable(clip_plane)

    glPopMatrix()


def draw_arrow(pos, d, radius, color, material, quality):
    pos2 = np.array(pos) + np.array(d)

    draw_cylinder(pos, pos2, radius, color, material, quality)

    ax, rx, ry = rotation_helper(d)

    glPushMatrix()
    glTranslatef(pos2[0], pos2[1], pos2[2])
    glRotatef(ax, rx, ry, 0.0)
    glutSolidCone(radius * 3, radius * 3, quality, quality)
    glPopMatrix()


# MOUSE EVENT MANAGER
class MouseFireEvent(object):

    """Event type of mouse button used for mouse callbacks.

    """

    ButtonPressed = 0
    FreeMotion = 1
    ButtonMotion = 2
    ButtonReleased = 3


class MouseButtonEvent(object):

    """Mouse event used for mouse callbacks. Stores button and callback.

    """

    def __init__(self, button, fireEvent, callback, positional=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.positional = positional


class MouseManager(object):

    """Handles mouse callbacks.

    """

    def __init__(self):
        self.mousePos = np.array([0, 0])
        self.mousePosOld = np.array([0, 0])
        self.mouseEventsPressed = []
        self.mouseEventsFreeMotion = []
        self.mouseEventsButtonMotion = []
        self.mouseEventsReleased = []
        self.mouseState = {}
        self.mouseState[GLUT_LEFT_BUTTON] = GLUT_UP
        self.mouseState[GLUT_MIDDLE_BUTTON] = GLUT_UP
        self.mouseState[GLUT_RIGHT_BUTTON] = GLUT_UP

    def register_button(self, mouseEvent):
        """Register mouse input callbacks.

        """
        if mouseEvent.fireEvent == MouseFireEvent.ButtonPressed:
            self.mouseEventsPressed.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonReleased:
            self.mouseEventsReleased.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.FreeMotion:
            self.mouseEventsFreeMotion.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonMotion:
            self.mouseEventsButtonMotion.append(mouseEvent)

    def mouse_click(self, button, state, x, y):
        self.mousePosOld = self.mousePos
        self.mousePos = np.array([x, y])

        self.mouseState[button] = state

        for me in self.mouseEventsPressed:
            if me.button == button and state == GLUT_DOWN:
                if me.positional:
                    me.callback(self.mousePos, self.mousePosOld)
                else:
                    me.callback()
        for me in self.mouseEventsReleased:
            if me.button == button and state == GLUT_UP:
                if me.positional:
                    me.callback(self.mousePos, self.mousePosOld)
                else:
                    me.callback()

    def mouse_move(self, x, y):
        self.mousePosOld = self.mousePos
        self.mousePos = np.array([x, y])

        for me in self.mouseEventsFreeMotion:
            me.callback(self.mousePos, self.mousePosOld, self.mouseState)

# KEYBOARD EVENT MANAGER


class KeyboardFireEvent(object):

    """Event type of button used for keyboard callbacks.

    """

    Pressed = 0
    Hold = 1
    Released = 2


class KeyboardButtonEvent(object):

    """Keyboard event used for keyboard callbacks. Stores button, event type and callback.

    """

    def __init__(self, button, fireEvent, callback, internal=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.internal = internal


class KeyboardManager(object):

    """Handles keyboard callbacks.

    """

    def __init__(self):
        self.pressedKeys = set([])
        self.keyStateOld = {}
        self.keyState = {}
        self.buttonEventsPressed = []
        self.buttonEventsHold = []
        self.buttonEventsReleased = []
        self.userCallbackStack = []

    def register_button(self, buttonEvent):
        """Register keyboard input callbacks.

        """
        if buttonEvent.fireEvent == KeyboardFireEvent.Pressed:
            self.buttonEventsPressed.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Hold:
            self.buttonEventsHold.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Released:
            self.buttonEventsReleased.append(buttonEvent)

    def callback_on_button(self, be, b):
        if be.button == b:
            if be.internal:
                be.callback()
            else:
                self.userCallbackStack.append(be.callback)

    def handle_input(self):
        removeKeys = set([])
        for b in self.pressedKeys:
            if self.keyStateOld[b] == 0 and self.keyState[b] == 1:
                for be in self.buttonEventsPressed:
                    self.callback_on_button(be, b)
                for be in self.buttonEventsHold:
                    self.callback_on_button(be, b)

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 1:
                for be in self.buttonEventsHold:
                    self.callback_on_button(be, b)

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 0:
                for be in self.buttonEventsReleased:
                    self.callback_on_button(be, b)
                removeKeys.add(b)

            self.keyStateOld[b] = self.keyState[b]

        self.pressedKeys = self.pressedKeys.difference(removeKeys)

    def keyboard_up(self, button):
        self.keyState[button] = 0  # Key up

    def keyboard_down(self, button):
        self.pressedKeys.add(button)
        self.keyState[button] = 1  # Key down
        if not button in self.keyStateOld.keys():
            self.keyStateOld[button] = 0

# CAMERA


class Camera(object):

    def __init__(self, camPos=np.array([0, 0, 1]), camTarget=np.array([0, 0, 0]), camRight=np.array([1.0, 0.0, 0.0]), moveSpeed=0.5, rotSpeed=0.001, globalRotSpeed=3.0, center=np.array([0, 0, 0]), updateLights=None):
        self.moveSpeed = moveSpeed
        self.lookSpeed = rotSpeed
        self.globalRotSpeed = globalRotSpeed

        self.center = center
        self.updateLights = updateLights

        self.modelview = np.identity(4, np.float32)

        t = camPos - camTarget
        r = np.linalg.norm(t)

        self.state_target = -t / r
        self.state_right = camRight / np.linalg.norm(camRight)
        self.state_up = np.cross(self.state_right, self.state_target)

        self.state_pos = np.array([0, 0, r])

        self.update_modelview()

    def move_forward(self):
        self.state_pos[2] += self.moveSpeed
        self.update_modelview()

    def move_backward(self):
        self.state_pos[2] -= self.moveSpeed
        self.update_modelview()

    def move_up(self):
        self.state_pos[1] += self.moveSpeed
        self.update_modelview()

    def move_down(self):
        self.state_pos[1] -= self.moveSpeed
        self.update_modelview()

    def move_left(self):
        self.state_pos[0] -= self.moveSpeed
        self.update_modelview()

    def move_right(self):
        self.state_pos[0] += self.moveSpeed
        self.update_modelview()

    def rotate_system_XL(self):
        self.rotate_camera_H(0.01 * self.globalRotSpeed)

    def rotate_system_XR(self):
        self.rotate_camera_H(-0.01 * self.globalRotSpeed)

    def rotate_system_YL(self):
        self.rotate_camera_R(0.01 * self.globalRotSpeed)

    def rotate_system_YR(self):
        self.rotate_camera_R(-0.01 * self.globalRotSpeed)

    def rotate_system_ZL(self):
        self.rotate_camera_V(0.01 * self.globalRotSpeed)

    def rotate_system_ZR(self):
        self.rotate_camera_V(-0.01 * self.globalRotSpeed)

    def rotate_camera(self, mousePos, mousePosOld, mouseButtonState):
        dm = mousePos - mousePosOld

        if mouseButtonState[GLUT_LEFT_BUTTON] == GLUT_DOWN:
            if dm[0] != 0:
                self.rotate_camera_H(dm[0] * 0.001 * self.globalRotSpeed)
            if dm[1] != 0:
                self.rotate_camera_V(dm[1] * 0.001 * self.globalRotSpeed)
        elif mouseButtonState[GLUT_RIGHT_BUTTON] == GLUT_DOWN:
            self.state_pos[0] -= 0.05 * dm[0] * self.moveSpeed
            self.state_pos[1] += 0.05 * dm[1] * self.moveSpeed
            self.update_modelview()
        elif mouseButtonState[GLUT_MIDDLE_BUTTON] == GLUT_DOWN:
            self.state_pos[2] += 0.05 * dm[1] * self.moveSpeed
            self.rotate_camera_R(dm[0] * 0.001 * self.globalRotSpeed)

    def normalize(self, vec):
        vec = self.normalized(vec)

    def normalized(self, vec):
        return vec / np.linalg.norm(vec)

    def get_camera_rotation_matrix(self, target_vec, up_vec):
        n = self.normalized(target_vec)
        u = self.normalized(up_vec)
        u = np.cross(u, target_vec)
        v = np.cross(n, u)
        m = np.identity(4, np.float32)
        m[0][0] = u[0]
        m[0][1] = v[0]
        m[0][2] = n[0]
        m[1][0] = u[1]
        m[1][1] = v[1]
        m[1][2] = n[1]
        m[2][0] = u[2]
        m[2][1] = v[2]
        m[2][2] = n[2]
        return m

    def rotate_vector(self, vec, ang, axe):
        sinhalf = sin(ang / 2)
        coshalf = cos(ang / 2)

        rx = axe[0] * sinhalf
        ry = axe[1] * sinhalf
        rz = axe[2] * sinhalf
        rw = coshalf

        rot = Quaternion(rx, ry, rz, rw)
        conq = rot.conjugated()
        w1 = conq.mult_v(vec)
        w = w1.mult_q(rot)

        vec[0] = w[0]
        vec[1] = w[1]
        vec[2] = w[2]

    def rotate_camera_R(self, da):
        self.rotate_vector(self.state_right, da, self.state_target)
        self.rotate_vector(self.state_up, da, self.state_target)
        self.update_modelview()

    def rotate_camera_V(self, da):
        self.rotate_vector(self.state_target, da, self.state_right)
        self.state_up = np.cross(self.state_right, self.state_target)
        self.update_modelview()

    def rotate_camera_H(self, da):
        self.rotate_vector(self.state_target, da, self.state_up)
        self.state_right = np.cross(self.state_target, self.state_up)
        self.update_modelview()

    def update_modelview(self):

        self.state_up /= np.linalg.norm(self.state_up)
        self.state_right /= np.linalg.norm(self.state_right)
        self.state_target /= np.linalg.norm(self.state_target)

        # Center Box
        trans = np.identity(4, np.float32)
        trans[3][0] = -self.center[0]
        trans[3][1] = -self.center[1]
        trans[3][2] = -self.center[2]

        # Camera rotation
        rotate_cam = self.get_camera_rotation_matrix(
            -self.state_target, self.state_up)

        # System translation
        trans_cam = np.identity(4, np.float32)
        trans_cam[3][0] = -self.state_pos[0]
        trans_cam[3][1] = -self.state_pos[1]
        trans_cam[3][2] = -self.state_pos[2]

        self.modelview = trans.dot(rotate_cam.dot(trans_cam))

        cXYZ = -1 * np.mat(self.modelview[:3, :3]) * \
            np.mat(self.modelview[3, :3]).T
        self.camPos = np.array([cXYZ[0, 0], cXYZ[1, 0], cXYZ[2, 0]])

        self.updateLights()


class Quaternion:

    def __init__(self, x, y, z, w):
        self.array = np.array([x, y, z, w], np.float32)

    def __getitem__(self, x):
        return self.array[x]

    def conjugated(self):
        return Quaternion(-self[0], -self[1], -self[2], self[3])

    def mult_v(self, v):
        w = - (self[0] * v[0]) - (self[1] * v[1]) - (self[2] * v[2])
        x = (self[3] * v[0]) + (self[1] * v[2]) - (self[2] * v[1])
        y = (self[3] * v[1]) + (self[2] * v[0]) - (self[0] * v[2])
        z = (self[3] * v[2]) + (self[0] * v[1]) - (self[1] * v[0])
        return Quaternion(x, y, z, w)

    def mult_q(self, q):
        w = - (self[3] * q[3]) - (self[0] * q[0]) - \
            (self[1] * q[1]) - (self[2] * q[2])
        x =   (self[0] * q[3]) + (self[3] * q[0]) + \
            (self[1] * q[2]) - (self[2] * q[1])
        y = (self[1] * q[3]) + (self[3] * q[1]) + (
            self[2] * q[0]) - (self[0] * q[2])
        z = (self[2] * q[3]) + (self[3] * q[2]) + (
            self[0] * q[1]) - (self[1] * q[0])
        return Quaternion(x, y, z, w)
