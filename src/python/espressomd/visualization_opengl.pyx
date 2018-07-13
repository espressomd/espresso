from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from math import *
import numpy as np
import os
import time
import collections
import sys
from threading import Thread
from matplotlib.pyplot import imsave
include "myconfig.pxi"
from copy import deepcopy
import espressomd
from espressomd.particle_data import ParticleHandle


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
                Draw xyz system axes.
    draw_nodes : :obj:`bool`, optional
                Draw node boxes.
    draw_cells : :obj:`bool`, optional
                Draw cell boxes.
    quality_particles : :obj:`int`, optional
                        The number of subdivisions for particle spheres.
    quality_bonds : :obj:`int`, optional
                    The number of subdivisions for cylindrical bonds.
    quality_arrows : :obj:`int`, optional
                     The number of subdivisions for external force arrows.
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
                        node: Color according to the node the particle is on.
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
    ext_force_arrows_type_scale : array_like :obj:`float`, optional
                                  List of scale factors of external force arrows for different particle types.
    ext_force_arrows_type_colors : array_like :obj:`float`, optional
                                   Colors of ext_force arrows for different particle types.
    ext_force_arrows_type_materials : array_like :obj:`float`, optional
                                      Materils of ext_force arrows for different particle types.
    ext_force_arrows_type_radii : array_like :obj:`float`, optional
                                   List of arrow radii for different particle types.
    force_arrows : :obj:`bool`, optional
                   Enables particle force visualization.
    force_arrows_type_scale : array_like :obj:`float`, optional
                              List of scale factors of particle force arrows for different particle types.
    force_arrows_type_colors : array_like :obj:`float`, optional
                               Colors of particle force arrows for different particle types.
    force_arrows_type_materials : array_like :obj:`float`, optional
                                  Materials of particle force arrows for different particle types.
    force_arrows_type_radii : array_like :obj:`float`, optional
                               List of arrow radii for different particle types.
    velocity_arrows : :obj:`bool`, optional
                       Enables particle velocity visualization.
    velocity_arrows_type_scale : array_like :obj:`float`, optional
                                 List of scale factors of particle velocity arrows for different particle types.
    velocity_arrows_type_colors : array_like :obj:`float`, optional
                                  Colors of particle velocity arrows for different particle types.
    velocity_arrows_type_materials : array_like :obj:`float`, optional
                                     Materials of particle velocity arrows for different particle types.
    velocity_arrows_type_radii : array_like :obj:`float`, optional
                                  List of arrow radii for different particle types.
    director_arrows : :obj:`bool`, optional
                       Enables particle director visualization.
    director_arrows_type_scale : :obj:`float`, optional
                             Scale factor of particle director arrows for different particle types.
    director_arrows_type_colors : array_like :obj:`float`, optional
                                  Colors of particle director arrows for different particle types.
    director_arrows_type_materials : array_like :obj:`float`, optional
                                     Materials of particle director arrows for different particle types.
    director_arrows_type_radii : array_like :obj:`float`, optional
                                  List of arrow radii for different particle types.
    drag_enabled : :obj:`bool`, optional
                   Enables mouse-controlled particles dragging (Default: False)
    drag_force : :obj:`bool`, optional
                 Factor for particle dragging

    LB_draw_nodes : :obj:`bool`, optional
                   Draws a lattice representation of the LB nodes that are no boundaries.
    LB_draw_node_boundaries : :obj:`bool`, optional
                             Draws a lattice representation of the LB nodes that are boundaries.
    LB_draw_boundaries : :obj:`bool`, optional
                        Draws the LB shapes.
    LB_draw_velocity_plane : :obj:`bool`, optional
                             Draws LB node velocity arrows specified by LB_plane_axis, LB_plane_dist, LB_plane_ngrid.
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
            'bright':       [0.9, 1.0, 0.8, 0.4, 1.0],
            'medium':       [0.6, 0.8, 0.2, 0.4, 1.0],
            'dark':         [0.4, 0.5, 0.1, 0.4, 1.0],
            'transparent1': [0.6, 0.8, 0.2, 0.5, 0.8],
            'transparent2': [0.6, 0.8, 0.2, 0.5, 0.4],
            'transparent3': [0.6, 0.8, 0.2, 0.5, 0.2],
            'rubber':   	[0, 0.4, 0.7, 0.078125, 1.0],
            'chrome':       [0.25, 0.4, 0.774597, 0.6, 1.0],
            'plastic': 		[0, 0.55, 0.7, 0.25, 1.0],
            'steel':		[0.25, 0.38, 0, 0.32, 1.0]
        }

        # DEFAULT PROPERTIES
        self.specs = {
            'window_size': [800, 800],
            'name': 'Espresso Visualization',

            'background_color': [0, 0, 0],

            'draw_fps': False,
            'draw_box': True,
            'draw_axis': True,
            'draw_nodes': False,
            'draw_cells': False,

            'update_fps': 30,

            'periodic_images': [0, 0, 0],

            'quality_particles': 20,
            'quality_bonds': 16,
            'quality_arrows': 16,
            'quality_constraints': 32,

            'close_cut_distance': 0.1,
            'far_cut_distance': 5,

            'camera_position': 'auto',
            'camera_target': 'auto',
            'camera_right': [1.0, 0.0, 0.0],

            'particle_coloring': 'auto',
            'particle_sizes': 'auto',
            'particle_type_colors': [[1, 1, 0], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 1], [1, 0.5, 0], [0.5, 0, 1]],
            'particle_type_materials': ['medium'],
            'particle_charge_colors': [[1, 0, 0], [0, 1, 0]],

            'draw_constraints': True,
            'rasterize_pointsize': 10,
            'rasterize_resolution': 75.0,
            'constraint_type_colors': [[0.5, 0.5, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0], [0, 0, 0.5], [0.5, 0, 0]],
            'constraint_type_materials': ['transparent1'],

            'draw_bonds': True,
            'bond_type_radius': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            'bond_type_colors': [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'bond_type_materials': ['medium'],

            'ext_force_arrows': False,
            'ext_force_arrows_type_scale': [1.0],
            'ext_force_arrows_type_colors': [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'ext_force_arrows_type_materials': ['transparent2'],
            'ext_force_arrows_type_radii': [0.2],

            'velocity_arrows': False,
            'velocity_arrows_type_scale': [1.0],
            'velocity_arrows_type_colors': [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'velocity_arrows_type_materials': ['transparent2'],
            'velocity_arrows_type_radii': [0.2],

            'force_arrows': False,
            'force_arrows_type_scale': [1.0],
            'force_arrows_type_colors': [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'force_arrows_type_materials': ['transparent2'],
            'force_arrows_type_radii': [0.2],

            'director_arrows': False,
            'director_arrows_type_scale': [1.0],
            'director_arrows_type_colors': [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1], [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'director_arrows_type_materials': ['transparent2'],
            'director_arrows_type_radii': [0.2],

            'LB_draw_nodes': False,
            'LB_draw_node_boundaries': False,
            'LB_draw_boundaries': False,

            'LB_draw_velocity_plane': False,
            'LB_plane_axis': 2,
            'LB_plane_dist': 0,
            'LB_plane_ngrid': 5,
            'LB_vel_scale': 1.0,
            'LB_arrow_color': [1, 1, 1],
            'LB_arrow_material': 'transparent1',
            'LB_arrow_quality': 16,

            'light_pos': 'auto',
            'light_colors': [[0.1, 0.1, 0.2], [0.9, 0.9, 0.9], [1.0, 1.0, 1.0]],
            'light_brightness': 1.0,
            'light_size': 'auto',

            'spotlight_enabled': False,
            'spotlight_colors': [[0.2, 0.2, 0.3], [0.5, 0.5, 0.5], [1.0, 1.0, 1.0]],
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
            self.specs['ext_force_arrows'] = False

        IF not ROTATION:
            self.specs['director_arrows'] = False

        IF not LB and not LB_GPU:
            self.specs['LB_draw_velocity_plane'] = False
            self.specs['LB_draw_boundaries'] = False
            self.specs['LB_draw_nodes'] = False
            self.specs['LB_draw_node_boundaries'] = False

        IF not LB_BOUNDARIES and not LB_BOUNDARIES_GPU:
            self.specs['LB_draw_boundaries'] = False
            self.specs['LB_draw_node_boundaries'] = False

        # ESPRESSO RELATED INITS THAT ARE KNOWN ONLY WHEN RUNNING THE
        # INTEGRATION LOOP ARE CALLED ONCE IN UPDATE LOOP:
        # CONSTRAINTS, NODE BOXES, CELL BOXES, CHARGE RANGE, BONDS

        # ESPRESSO RELATED INITS THAT ARE KNOWN ON INSTANTIATION GO HERE:

        # CONTENT OF PARTICLE DATA
        self.has_particle_data = {}
        self.has_particle_data['velocity'] = self.specs['velocity_arrows']
        self.has_particle_data['force'] = self.specs['force_arrows']
        self.has_particle_data['ext_force'] = self.specs['ext_force_arrows'] or self.specs['drag_enabled']
        IF ELECTROSTATICS:
            self.has_particle_data['charge'] = self.specs['particle_coloring'] == 'auto' or self.specs['particle_coloring'] == 'charge'
        ELSE:
            self.has_particle_data['charge'] = False
        self.has_particle_data['director'] = self.specs['director_arrows']
        self.has_particle_data['node'] = self.specs['particle_coloring'] == 'node'

        # PARTICLE INFO OF HIGHLIGHTED PARTICLE: COLLECT PARTICLE ATTRIBUTES
        self.highlighted_particle = {}
        self.particle_attributes = []
        for d in dir(ParticleHandle):
            if type(getattr(ParticleHandle, d)) == type(ParticleHandle.pos):
                if not d in ["pos_folded"]:
                    self.particle_attributes.append(d)
        self.max_len_attr = max([len(a) for a in self.particle_attributes])

        # FIXED COLORS FROM INVERSE BACKGROUND COLOR FOR GOOD CONTRAST
        self.invBackgroundCol = np.array([1 - self.specs['background_color'][0], 1 -
                                          self.specs['background_color'][1], 1 - self.specs['background_color'][2], 1.0])

        self.node_box_color = np.copy(self.invBackgroundCol)
        self.node_box_color[0] += 0.5 * (0.5 - self.node_box_color[0])

        self.cell_box_color = np.copy(self.invBackgroundCol)
        self.cell_box_color[1] += 0.5 * (0.5 - self.cell_box_color[1])

        self.lb_box_color = np.copy(self.invBackgroundCol)
        self.lb_box_color[2] = 0.5

        self.lb_box_color_boundary = np.copy(self.invBackgroundCol)
        self.lb_box_color_boundary[1] = 0.5

        self.text_color = np.copy(self.invBackgroundCol)

        # INCREASE LINE THICKNESS IF NODE/CELL BOX IS ENABLED
        self.line_width_fac = 1.0
        if self.specs['draw_nodes']:
            self.line_width_fac += 0.5
        if self.specs['draw_cells']:
            self.line_width_fac += 0.5

        # HAS PERIODIC IMAGES
        self.has_images = any(i != 0 for i in self.specs['periodic_images'])

        # INITS
        self.system = system

        self.last_T = -1
        self.last_box_l = self.system.box_l
        self.fps_last = 0
        self.fps = 0
        self.fps_count = 0

        self.glutMainLoop_started = False
        self.screenshot_initialized = False
        self.hasParticleData = False
        self.quit_savely = False
        self.paused = False
        self.take_screenshot = False
        self.screenshot_captured = False

        self.keyboardManager = KeyboardManager()
        self.mouseManager = MouseManager()
        self.camera = Camera()
        self._init_camera()

        self.timers = []
        self.particles = {}

        self.update_elapsed = 0
        self.update_timer = 0
        self.draw_elapsed = 0
        self.draw_timer = 0

        # LIST OF [[px,py],[string]] FOR USER DEFINED TEXT
        self.user_texts = []

    def update_system_info(self):

        # SYSTEM INFORMATION
        self.system_info = {}
        self.system_info['Actors'] = []
        self.system_info['Non-bonded interactions'] = []
        self.system_info['Bonded interactions'] = [
            b for b in self.system.bonded_inter]
        self.system_info['Constraints'] = []
        self.system_info['Thermostat'] = self.system.thermostat.get_state()

        if len(self.system_info['Bonded interactions']) == 0:
            self.system_info['Bonded interactions'].append('None')

        # ACTORS
        for a in self.system.actors:
            self.system_info['Actors'].append(str(a))
        if len(self.system_info['Actors']) == 0:
            self.system_info['Actors'].append('None')

        # COLLECT ALL TYPES FROM PARTICLES AND CONSTRAINTS
        all_types = set(self.system.part[:].type)
        constraint_types = []
        for c in self.system.constraints[:]:
            constraint_types.append(c.get_parameter('particle_type'))
        all_types.update(constraint_types)

        # COLLECT ALL ACTIVCE NONBONDED INTERACTIONS
        all_non_bonded_inters = [x for x in dir(self.system.non_bonded_inter[0, 0]) if not x.startswith(
            '__') and not x == 'type1' and not x == 'type2']
        for t1 in all_types:
            for t2 in all_types:
                for check_nb in all_non_bonded_inters:
                    nb = getattr(
                        self.system.non_bonded_inter[t1, t2], check_nb)
                    if not nb == None and nb.is_active():
                        self.system_info['Non-bonded interactions'].append(
                            [t1, t2, nb.type_name(), nb.get_params()])
        if len(self.system_info['Non-bonded interactions']) == 0:
            self.system_info['Non-bonded interactions'].append('None')

        # COLLECT CONSTRAINTS
        for c in self.system.constraints[:]:
            co = c.get_params()
            co_s = c.get_parameter('shape')
            co['shape'] = [co_s.name(), co_s.get_params()]
            self.system_info['Constraints'].append(co)
        if len(self.system_info['Constraints']) == 0:
            self.system_info['Constraints'].append('None')

    def register_callback(self, cb, interval=1000):
        """Register timed callbacks.
        """

        self.timers.append((int(interval), cb))

    def screenshot(self, path):
        """Renders the current state and into an image file at path with dimensions of
        specs['window_size'].  """

        # ON FIRST CALL: INIT AND CREATE BUFFERS
        if not self.screenshot_initialized:
            self.screenshot_initialized = True
            self._init_opengl()

            # CREATE BUFFERS THAT CAN BE LARGER THAN THE SCREEN
            # FRAME BUFFER
            fbo = glGenFramebuffers(1)
            glBindFramebuffer(GL_FRAMEBUFFER, fbo)
            # COLOR BUFFER
            rbo = glGenRenderbuffers(1)
            glBindRenderbuffer(GL_RENDERBUFFER, rbo)
            glRenderbufferStorage(
                GL_RENDERBUFFER, GL_RGB, self.specs['window_size'][0], self.specs['window_size'][1])
            glFramebufferRenderbuffer(
                GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rbo)
            # DEPTH BUFFER
            dbo = glGenRenderbuffers(1)
            glBindRenderbuffer(GL_RENDERBUFFER, dbo)
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT,
                                  self.specs['window_size'][0], self.specs['window_size'][1])
            glFramebufferRenderbuffer(
                GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, dbo)

            self._reshape_window(
                self.specs['window_size'][0], self.specs['window_size'][1])
            glutHideWindow()

        # INIT AND UPDATE ESPRESSO
        self._init_espresso_visualization()
        self._initial_espresso_updates()

        # DRAW
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadMatrixf(self.camera.modelview)
        self._draw_system()

        # READ THE PIXES
        glReadBuffer(GL_COLOR_ATTACHMENT0)
        data = glReadPixels(
            0, 0, self.specs['window_size'][0], self.specs['window_size'][1], GL_RGB, GL_FLOAT)

        # RESHAPE THE DATA
        data = np.flipud(data.reshape((data.shape[1], data.shape[0], 3)))

        # SAVE TO IMAGE
        imsave(path, data)

    def run(self, integ_steps=1):
        """Convenience method with a simple integration thread.
        """

        def main():
            while True:

                self.update()

                if self.paused:
                    # sleep(0) is worse
                    time.sleep(0.0001)
                else:
                    try:
                        self.system.integrator.run(integ_steps)
                    except Exception as e:
                        print(e)
                        os._exit(1)

        t = Thread(target=main)
        t.daemon = True
        t.start()

        self.start()

    def start(self):
        """The blocking start method.
        """

        self._init_opengl()
        self._init_espresso_visualization()
        self._init_controls()
        self._init_openGL_callbacks()
        self._init_timers()
        self._init_camera()

        # START THE BLOCKING MAIN LOOP
        self.glutMainLoop_started = True
        glutMainLoop()

    def update(self):
        """Update method to be called after integration.
        Changes of espresso system can only happen here.

        """
        if self.glutMainLoop_started:

            # UPDATE ON STARTUP
            if not self.hasParticleData:
                self._initial_espresso_updates()
                self.hasParticleData = True

            # UPDATES
            self.update_elapsed += (time.time() - self.update_timer)
            if self.update_elapsed > 1.0 / self.specs['update_fps']:
                self.update_elapsed = 0

                # ES UPDATES WHEN SYSTEM HAS PROPAGATED. ALSO UPDATE ON PAUSE FOR PARTICLE INFO
                if self.paused or not self.last_T == self.system.time:
                    self.last_T = self.system.time
                    self._update_particles()

                    # LB UPDATE
                    if self.specs['LB_draw_velocity_plane']:
                        self._update_lb_velocity_plane()

                    # BOX_L CHANGED
                    if not (self.last_box_l == self.system.box_l).all():
                        self.last_box_l = np.copy(self.system.box_l)
                        self._box_size_dependence()

                # KEYBOARD CALLBACKS MAY CHANGE ESPRESSO SYSTEM PROPERTIES,
                # ONLY SAVE TO CHANGE HERE
                for c in self.keyboardManager.userCallbackStack:
                    c()
                self.keyboardManager.userCallbackStack = []

            self.update_timer = time.time()

            # DRAG PARTICLES
            if self.specs['drag_enabled']:
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

    # FETCH DATA ON STARTUP
    def _initial_espresso_updates(self):
        self._update_particles()
        if self.has_particle_data['charge']:
            self._update_charge_color_range()
        self._update_bonds()
        if self.specs['draw_constraints']:
            self._update_constraints()
        if self.specs['draw_cells'] or self.specs['draw_nodes']:
            self._update_nodes()
        if self.specs['draw_cells']:
            self._update_cells()

    # GET THE PARTICLE DATA
    def _update_particles(self):

        self.particles['pos'] = self.system.part[:].pos_folded
        self.particles['type'] = self.system.part[:].type

        if self.has_particle_data['velocity']:
            self.particles['velocity'] = self.system.part[:].v

        if self.has_particle_data['force']:
            self.particles['force'] = self.system.part[:].f

        if self.has_particle_data['ext_force']:
            self.particles['ext_force'] = self.system.part[:].ext_force

        if self.has_particle_data['charge']:
            self.particles['charge'] = self.system.part[:].q

        if self.has_particle_data['director']:
            self.particles['director'] = self.system.part[:].director

        if self.has_particle_data['node']:
            self.particles['node'] = self.system.part[:].node

        if self.infoId != -1:
            for attr in self.particle_attributes:
                self.highlighted_particle[attr] = getattr(
                    self.system.part[self.infoId], attr)

    def _update_lb_velocity_plane(self):
        if self.lb_is_cpu:
            self._update_lb_velocity_plane_CPU()
        else:
            self._update_lb_velocity_plane_GPU()

    def _update_lb_velocity_plane_CPU(self):
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

    def _update_lb_velocity_plane_GPU(self):
        agrid = self.lb_params['agrid']
        ng = self.specs['LB_plane_ngrid']
        col_pos = []
        for xi in xrange(ng):
            for xj in xrange(ng):
                p = np.array((self.lb_plane_p + xi * 1.0 / ng * self.lb_plane_b1 +
                              xj * 1.0 / ng * self.lb_plane_b2) % self.system.box_l)
                col_pos.append(p)

        lb_vels = self.lb.get_interpolated_fluid_velocity_at_positions(
            np.array(col_pos))
        self.lb_plane_vel = []
        for p, v in zip(col_pos, lb_vels):
            self.lb_plane_vel.append([p, v])

    def _edges_from_pn(self, p, n, diag):
        v1, v2 = self._get_tangents(n)

        edges = []
        edges.append(p + diag * v1)
        edges.append(p + diag * v2)
        edges.append(p - diag * v1)
        edges.append(p - diag * v2)
        return edges

    def _update_cells(self):
        self.cell_box_origins = []
        cell_system_state = self.system.cell_system.get_state()
        self.cell_size = cell_system_state['cell_size']
        for i in range(cell_system_state['cell_grid'][0]):
            for j in range(cell_system_state['cell_grid'][1]):
                for k in range(cell_system_state['cell_grid'][2]):
                    self.cell_box_origins.append(
                        np.array([i, j, k]) * self.cell_size)

    def _update_nodes(self):
        self.node_box_origins = []
        self.local_box_l = self.system.cell_system.get_state()['local_box_l']
        for i in range(self.system.cell_system.node_grid[0]):
            for j in range(self.system.cell_system.node_grid[1]):
                for k in range(self.system.cell_system.node_grid[2]):
                    self.node_box_origins.append(
                        np.array([i, j, k]) * self.local_box_l)

    # GET THE _update_constraints DATA
    def _update_constraints(self):

        box_diag = pow(pow(self.system.box_l[0], 2) + pow(
            self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)

        self.shapes = collections.defaultdict(list)

        # Collect shapes and interaction type (for coloring) from constraints
        primitive_shapes = ['Shapes::Wall', 'Shapes::Cylinder', 'Shapes::Ellipsoid',
                            'Shapes::SimplePore', 'Shapes::Sphere', 'Shapes::SpheroCylinder']

        coll_shape_obj = collections.defaultdict(list)
        for c in self.system.constraints:
            if type(c) == espressomd.constraints.ShapeBasedConstraint:
                t = c.get_parameter('particle_type')
                s = c.get_parameter('shape')
                n = s.name()
                if n in primitive_shapes:
                    coll_shape_obj[n].append([s, t])
                else:
                    coll_shape_obj['Shapes::Misc'].append([s, t])

        if self.specs['LB_draw_boundaries']:
            ni = 0
            for c in self.system.lbboundaries:
                if type(c) == espressomd.ekboundaries.EKBoundary:
                    t = ni
                    ni += 1
                    s = c.get_parameter('shape')
                    n = s.name()
                    if n in primitive_shapes:
                        coll_shape_obj[n].append([s, t])
                    else:
                        coll_shape_obj['Shapes::Misc'].append([s, t])

        for s in coll_shape_obj['Shapes::Wall']:
            d = s[0].get_parameter('dist')
            n = s[0].get_parameter('normal')
            edges = self._edges_from_pn(d * np.array(n), n, 2 * box_diag)
            self.shapes['Shapes::Wall'].append([edges, s[1]])

        for s in coll_shape_obj['Shapes::Cylinder']:
            pos = np.array(s[0].get_parameter('center'))
            a = np.array(s[0].get_parameter('axis'))
            l = s[0].get_parameter('length')
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::Cylinder'].append(
                [pos - a / np.linalg.norm(a) * l * 0.5, pos + a / np.linalg.norm(a) * l * 0.5, r, s[1]])

        for s in coll_shape_obj['Shapes::Ellipsoid']:
            pos = np.array(s[0].get_parameter('center'))
            a = np.array(s[0].get_parameter('a'))
            b = np.array(s[0].get_parameter('b'))
            c = np.array(s[0].get_parameter('b'))
            self.shapes['Shapes::Ellipsoid'].append([pos, a, b, c, s[1]])

        for s in coll_shape_obj['Shapes::Sphere']:
            pos = np.array(s[0].get_parameter('center'))
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::Sphere'].append([pos, r, s[1]])

        for s in coll_shape_obj['Shapes::SimplePore']:
            center = np.array(s[0].get_parameter('center'))
            axis = np.array(s[0].get_parameter('axis'))
            length = np.array(s[0].get_parameter('length'))
            radius = np.array(s[0].get_parameter('radius'))
            smoothing_radius = np.array(s[0].get_parameter('smoothing_radius'))
            self.shapes['Shapes::SimplePore'].append(
                [center, axis, length, radius, smoothing_radius, s[1]])

        for s in coll_shape_obj['Shapes::SpheroCylinder']:
            pos = np.array(s[0].get_parameter('center'))
            a = np.array(s[0].get_parameter('axis'))
            l = s[0].get_parameter('length')
            r = s[0].get_parameter('radius')
            self.shapes['Shapes::SpheroCylinder'].append(
                [pos - a / np.linalg.norm(a) * l * 0.5, pos + a / np.linalg.norm(a) * l * 0.5, r, s[1]])

        for s in coll_shape_obj['Shapes::Misc']:
            self.shapes['Shapes::Misc'].append(
                [self._rasterize_brute_force(s[0]), s[1]])

    def _get_tangents(self, n):
        n = np.array(n)
        v1 = np.random.randn(3)
        v1 -= v1.dot(n) * n / np.linalg.norm(n)**2
        v2 = np.cross(n, v1)
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        return v1, v2

    def _rasterize_brute_force(self, shape):
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
    def _update_bonds(self):
        if self.specs['draw_bonds']:
            self.bonds = []
            for i, p in enumerate(self.system.part):
                bs = p.bonds
                for b in bs:
                    t = b[0].type_number()
                    # b[0]: Bond, b[1:] Partners
                    for p in b[1:]:
                        self.bonds.append([i, p, t])

    def _draw_text(self, x,  y, text, color, font=GLUT_BITMAP_9_BY_15):
        glColor(color)
        glWindowPos2f(x, y)
        for ch in text:
            glutBitmapCharacter(font, ctypes.c_int(ord(ch)))

    # DRAW CALLED AUTOMATICALLY FROM GLUT DISPLAY FUNC
    def _draw_system(self):
        if self.specs['LB_draw_velocity_plane']:
            self._draw_lb_vel()

        if self.specs['draw_axis']:
            axis_fac = 0.2
            axis_r = np.min(self.system.box_l) / 50.0
            draw_arrow([0, 0, 0], [self.system.box_l[0] * axis_fac, 0, 0], axis_r, [
                       1, 0, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, self.system.box_l[2] * axis_fac, 0], axis_r, [
                       0, 1, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, 0, self.system.box_l[2] * axis_fac], axis_r, [
                       0, 0, 1, 1], self.materials['chrome'], self.specs['quality_arrows'])

        if self.specs['draw_bonds']:
            self._draw_bonds()
        self._draw_system_particles()

        if self.specs['draw_constraints']:
            self._draw_constraints()

        if self.specs['draw_box']:
            self._draw_system_box()
        if self.specs['draw_nodes']:
            self._draw_nodes()
        if self.specs['draw_cells']:
            self._draw_cells()
        if self.specs['LB_draw_nodes'] or self.specs['LB_draw_node_boundaries']:
            self._draw_lb_grid()

    def _draw_system_box(self):
        draw_box([0, 0, 0], self.system.box_l,
                 self.invBackgroundCol, self.materials['medium'], 2.0 * self.line_width_fac)

    def _draw_nodes(self):
        for n in self.node_box_origins:
            draw_box(n, self.local_box_l, self.node_box_color, self.materials['transparent1'],
                     1.5 * self.line_width_fac)

    def _draw_cells(self):
        for n in self.node_box_origins:
            for c in self.cell_box_origins:
                draw_box(c + n, self.cell_size, self.cell_box_color, self.materials['transparent1'],
                         0.75 * self.line_width_fac)

    def _draw_lb_grid(self):
        a = self.lb_params['agrid']
        cell_size = np.array([a] * 3)
        dims = np.rint(np.array(self.system.box_l) / a)
        for i in range(int(dims[0])):
            for j in range(int(dims[1])):
                for k in range(int(dims[2])):
                    n = np.array([i, j, k]) * cell_size
                    if self.specs['LB_draw_node_boundaries'] and self.lb[i, j, k].boundary:
                        draw_box(n, cell_size, self.lb_box_color_boundary,
                                 self.materials['transparent2'], 5.0)
                    if self.specs['LB_draw_nodes'] and not self.lb[i, j, k].boundary:
                        draw_box(n, cell_size, self.lb_box_color,
                                 self.materials['transparent2'], 1.5)

    def _draw_constraints(self):

        # CLIP BORDERS OF SIMULATION BOX
        for i in range(6):
            glEnable(GL_CLIP_PLANE0 + i)
            glClipPlane(GL_CLIP_PLANE0 + i, self.box_eqn[i])

        # NEEDS ADDITIONAL CLIP PLANES
        for s in self.shapes['Shapes::SimplePore']:
            draw_simple_pore(s[0], s[1], s[2], s[3], s[4], max(self.system.box_l), self._modulo_indexing(self.specs['constraint_type_colors'], s[5]),
                             self.materials[self._modulo_indexing(self.specs['constraint_type_materials'], s[5])], self.specs['quality_constraints'])

        # NEEDS ADDITIONAL CLIP PLANES
        for s in self.shapes['Shapes::SpheroCylinder']:
            draw_sphero_cylinder(
                s[0], s[1], s[2], self._modulo_indexing(
                    self.specs['constraint_type_colors'], s[3]),
                self.materials[self._modulo_indexing(self.specs['constraint_type_materials'], s[3])], self.specs['quality_constraints'])

        # RESET CLIP BORDERS
        for i in range(6):
            glEnable(GL_CLIP_PLANE0 + i)
            glClipPlane(GL_CLIP_PLANE0 + i, self.box_eqn[i])

        for s in self.shapes['Shapes::Ellipsoid']:
            draw_ellipsoid(s[0], s[1], s[2], s[3], self._modulo_indexing(self.specs['constraint_type_colors'], s[4]),
                           self.materials[self._modulo_indexing(self.specs['constraint_type_materials'], s[4])], self.specs['quality_constraints'])

        for s in self.shapes['Shapes::Sphere']:
            draw_sphere(s[0], s[1], self._modulo_indexing(self.specs['constraint_type_colors'], s[2]), self.materials[self._modulo_indexing(
                self.specs['constraint_type_materials'], s[2])], self.specs['quality_constraints'])

        for s in self.shapes['Shapes::Wall']:
            draw_plane(
                s[0], self._modulo_indexing(
                    self.specs['constraint_type_colors'], s[1]),
                self.materials[self._modulo_indexing(self.specs['constraint_type_materials'], s[1])])

        for s in self.shapes['Shapes::Cylinder']:
            draw_cylinder(s[0], s[1], s[2], self._modulo_indexing(self.specs['constraint_type_colors'], s[3]), self.materials[self._modulo_indexing(
                self.specs['constraint_type_materials'], s[3])], self.specs['quality_constraints'], True)

        for s in self.shapes['Shapes::Misc']:
            draw_points(s[0], self.specs['rasterize_pointsize'],  self._modulo_indexing(
                self.specs['constraint_type_colors'], s[1]), self.materials[self._modulo_indexing(self.specs['constraint_type_materials'], s[1])])

        for i in range(6):
            glDisable(GL_CLIP_PLANE0 + i)

    def _determine_radius(self, ptype):
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
                radius = self._modulo_indexing(
                    self.specs['particle_sizes'], ptype)
            except BaseException:
                radius = self.radiusByLJ(ptype)
        return radius

    def _draw_system_particles(self, colorById=False):
        pIds = range(len(self.particles['pos']))
        ptype = -1
        reset_material = False

        for pid in pIds:
            ptype_last = ptype
            ptype = int(self.particles['type'][pid])

            # Only change material if type/charge has changed, colorById or material was resetted by arrows
            if reset_material or colorById or not ptype == ptype_last or pid == self.dragId or pid == self.infoId or self.specs['particle_coloring'] == 'node':
                reset_material = False

                radius = self._determine_radius(ptype)

                m = self._modulo_indexing(
                    self.specs['particle_type_materials'], ptype)
                material = self.materials[m]

                if colorById:
                    color = self._id_to_fcolor(pid)
                    glColor(color)
                else:
                    if self.specs['particle_coloring'] == 'auto':
                        # Color auto: Charge then Type
                        if self.has_particle_data['charge'] and self.particles['charge'][pid] != 0:
                            color = self._color_by_charge(
                                self.particles['charge'][pid])
                            reset_material = True
                        else:
                            color = self._modulo_indexing(
                                self.specs['particle_type_colors'], ptype)
                    elif self.specs['particle_coloring'] == 'charge':
                        color = self._color_by_charge(
                            self.particles['charge'][pid])
                        reset_material = True
                    elif self.specs['particle_coloring'] == 'type':
                        color = self._modulo_indexing(
                            self.specs['particle_type_colors'], ptype)
                    elif self.specs['particle_coloring'] == 'node':
                        color = self._modulo_indexing(
                            self.specs['particle_type_colors'], self.particles['node'][pid])

                    # Invert color of highlighted particle
                    if pid == self.dragId or pid == self.infoId:
                        reset_material = True
                        color = [1 - color[0], 1 - color[1],
                                 1 - color[2]]

                    set_solid_material(color, material)

                # Create a new display list, used until next material/color change
                glNewList(self.dl_sphere, GL_COMPILE)
                glutSolidSphere(
                    radius, self.specs['quality_particles'], self.specs['quality_particles'])
                glEndList()

            self._redraw_sphere(
                self.particles['pos'][pid], radius, self.specs['quality_particles'])

            if self.has_images:
                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                self._redraw_sphere(
                                    self.particles['pos'][pid] + (imx * self.imPos[0] + imy * self.imPos[1] + imz * self.imPos[2]), radius, self.specs['quality_particles'])

            IF EXTERNAL_FORCES:
                if self.specs['ext_force_arrows'] or pid == self.dragId:
                    if any(v != 0 for v in self.particles['ext_force'][pid]):
                        if pid == self.dragId:
                            sc = 1
                        else:
                            sc = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_scale'], ptype)
                        if sc > 0:
                            arrow_col = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_colors'], ptype)
                            arrow_radius = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_radii'], ptype)
                            draw_arrow(self.particles['pos'][pid], np.array(
                                self.particles['ext_force'][pid]) * sc, arrow_radius, arrow_col, self.materials['chrome'], self.specs['quality_arrows'])
                            reset_material = True

            if self.specs['velocity_arrows']:
                self._draw_arrow_property(pid, ptype, self.specs['velocity_arrows_type_scale'], self.specs[
                    'velocity_arrows_type_colors'], self.specs['velocity_arrows_type_radii'], 'velocity')
                reset_material = True

            if self.specs['force_arrows']:
                self._draw_arrow_property(pid, ptype, self.specs['force_arrows_type_scale'],
                                          self.specs['force_arrows_type_colors'], self.specs['force_arrows_type_radii'], 'force')
                reset_material = True

            if self.specs['director_arrows']:
                self._draw_arrow_property(pid, ptype, self.specs['director_arrows_type_scale'], self.specs[
                    'director_arrows_type_colors'], self.specs['director_arrows_type_radii'], 'director')
                reset_material = True

    def _draw_arrow_property(self, pid, ptype, type_scale, type_colors, type_radii, prop):
        sc = self._modulo_indexing(type_scale, ptype)
        if sc > 0:
            v = self.particles[prop][pid]
            col = self._modulo_indexing(type_colors, ptype)
            radius = self._modulo_indexing(type_radii, ptype)
            draw_arrow(self.particles['pos'][pid], np.array(
                v, dtype=float) * sc, radius, col, self.materials['chrome'], self.specs['quality_arrows'])

    def _draw_bonds(self):
        pIds = range(len(self.particles['pos']))
        b2 = self.system.box_l[0] / 2.0
        box_l2_sqr = pow(b2, 2.0)
        for b in self.bonds:
            col = self._modulo_indexing(self.specs['bond_type_colors'], b[2])
            mat = self.materials[self._modulo_indexing(
                self.specs['bond_type_materials'], b[2])]
            radius = self._modulo_indexing(
                self.specs['bond_type_radius'], b[2])
            d = self.particles['pos'][b[0]] - self.particles['pos'][b[1]]
            bondLen_sqr = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]

            if bondLen_sqr < box_l2_sqr:
                # BOND COMPLETELY INSIDE BOX
                draw_cylinder(self.particles['pos'][b[0]], self.particles['pos'][b[1]], radius,
                              col, mat, self.specs['quality_bonds'])
                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                im = np.array([imx, imy, imz])
                                draw_cylinder(self.particles['pos'][b[0]] + im * self.imPos[dim], self.particles['pos'][b[1]] +
                                              im * self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])
            else:
                # SPLIT BOND
                l = self.particles['pos'][b[0]] - self.particles['pos'][b[1]]
                l0 = self.particles['pos'][b[0]]
                hits = 0
                for i in range(6):
                    lineBoxNDot = float(np.dot(l, self.box_n[i]))
                    if lineBoxNDot == 0:
                        continue
                    s = l0 - \
                        np.dot(l0 - self.box_p[i],
                               self.box_n[i]) / lineBoxNDot * l
                    if self._is_inside_box(s):
                        if lineBoxNDot < 0:
                            s0 = s
                        else:
                            s1 = s
                        hits += 1
                        if hits >= 2:
                            break
                if hits >= 2:
                    draw_cylinder(self.particles['pos'][b[0]], s0, radius, col,
                                  mat, self.specs['quality_bonds'])
                    draw_cylinder(self.particles['pos'][b[1]], s1, radius, col,
                                  mat, self.specs['quality_bonds'])

                    for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                        for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                            for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                                if imx != 0 or imy != 0 or imz != 0:
                                    im = np.array([imx, imy, imz])
                                    draw_cylinder(self.particles['pos'][b[0]] + im * self.imPos[dim], s0 + im *
                                                  self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])
                                    draw_cylinder(self.particles['pos'][b[1]] + im * self.imPos[dim], s1 + im *
                                                  self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])

    def _redraw_sphere(self, pos, radius, quality):
        glPushMatrix()
        glTranslatef(pos[0], pos[1], pos[2])
        glCallList(self.dl_sphere)
        glPopMatrix()

    # HELPER TO DRAW PERIODIC BONDS
    def _is_inside_box(self, p):
        eps = 1e-5
        for i in range(3):
            if p[i] < -eps or p[i] > eps + self.last_box_l[i]:
                return False
        return True

    # ARROWS IN A PLANE FOR LB VELOCITIES
    def _draw_lb_vel(self):

        for lbl in self.lb_plane_vel:
            p = lbl[0]
            v = lbl[1]
            c = np.linalg.norm(v)
            draw_arrow(p, v * self.specs['LB_vel_scale'], self.lb_arrow_radius, self.specs['LB_arrow_color'],
                       self.materials[self.specs['LB_arrow_material']], self.specs['LB_arrow_quality'])

    # USE MODULO IF THERE ARE MORE PARTICLE TYPES THAN TYPE DEFINITIONS FOR
    # COLORS, MATERIALS ETC..
    def _modulo_indexing(self, l, t):
        return l[t % len(l)]

    # FADE PARTICE CHARGE COLOR FROM WHITE (q=0) to PLUSCOLOR (q=q_max) RESP
    # MINUSCOLOR (q=q_min)
    def _color_by_charge(self, q):
        if q < 0:
            c = 1.0 * q / self.minq
            return np.array(self.specs['particle_charge_colors'][0]) * c + (1 - c) * np.array([1, 1, 1])
        else:
            c = 1.0 * q / self.maxq

            return np.array(self.specs['particle_charge_colors'][1]) * c + (1 - c) * np.array([1, 1, 1])

    # ON INITIALIZATION, CHECK q_max/q_min
    def _update_charge_color_range(self):
        if len(self.particles['charge'][:]) > 0:
            self.minq = min(self.particles['charge'][:])
            self.maxq = max(self.particles['charge'][:])

    def _handle_screenshot(self):
        if self.take_screenshot:
            self.take_screenshot = False
            data = glReadPixels(
                0, 0, self.specs['window_size'][0], self.specs['window_size'][1], GL_RGB, GL_FLOAT)
            scriptname = os.path.splitext(sys.argv[0])[0]

            i = 0
            while os.path.exists("{}_{num:04d}.png".format(scriptname, num=i)):
                i += 1
            fname = "{}_{num:04d}.png".format(scriptname, num=i)

            data = np.flipud(data.reshape((data.shape[1], data.shape[0], 3)))
            imsave(fname, data)

            self.screenshot_captured = True
            self.screenshot_capture_time = time.time()
            self.screenshot_capture_txt = "Saved screenshot {}".format(fname)

    def _display_all(self):

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glLoadMatrixf(self.camera.modelview)

        self._set_camera_spotlight()

        self._draw_system()
        self._draw_texts()

        glutSwapBuffers()

        self._handle_screenshot()

    def _draw_texts(self):

        # DRAW USER TEXT
        for ut in self.user_texts:
            p = ut[0]
            t = ut[1]
            self._draw_text(p[0], p[1], t, self.text_color)

        # DRAW FPS TEXT
        if self.specs['draw_fps']:
            t = time.time()
            if t - self.fps_last > 1.0:
                self.fps_last = t
                self.fps = self.fps_count
                self.fps_count = 0

            self._draw_text(10, 10, "{} fps".format(
                self.fps), self.text_color)
            self._draw_text(
                10, 30, "{} ms/frame".format(1000.0 / self.fps),  self.text_color)
            self.fps_count += 1

        # DRAW PARTICLE INFO
        if self.show_system_info:
            self._draw_sysinfo_dict(self.system_info)
        elif self.infoId != -1:
            self._draw_particle_dict(
                self.highlighted_particle, self.max_len_attr)

        # INDICATE SCREEN CAPTURE
        if self.screenshot_captured and not self.take_screenshot:
            col = np.array(self.text_color)
            ts = time.time() - self.screenshot_capture_time
            fadetime = 2.0
            col[3] = 1.0 - ts / fadetime
            if ts > fadetime:
                self.screenshot_captured = False
            else:
                self._draw_text(self.specs['window_size'][0] - len(self.screenshot_capture_txt) * 9.0 - 15,
                                self.specs['window_size'][1] - 15, self.screenshot_capture_txt, col)

    def _draw_sysinfo_dict(self, d):
        y = 0
        for k, v in d.items():
            # CATEGORY TITLE
            self._draw_text(
                10, self.specs['window_size'][1] - 10 - 15 * y, str(k) + ":", self.text_color)
            # ITEM LIST
            for item in v:
                txt = str(item)
                # NUMBER OF LINES
                nl = int(len(txt) * 9.0 / self.specs['window_size'][0]) + 1
                ch_per_line = int(self.specs['window_size'][0] / 9) - 4
                if ch_per_line < 20:
                    break
                ls = 0
                for li in range(nl):
                    y += 1
                    ltxt = txt[ls:ls + ch_per_line]
                    self._draw_text(
                        30, self.specs['window_size'][1] - 10 - 15 * y, ltxt, self.text_color)
                    ls += ch_per_line
            y += 1.5

    def _draw_particle_dict(self, d, maxlen):
        y = 0
        for k, v in d.items():
            txt = "{} {} {}".format(
                k, (maxlen - len(k)) * ' ',  v)
            self._draw_text(
                10, self.specs['window_size'][1] - 10 - 15 * y, txt, self.text_color)
            y += 1

    # CALLED ION WINDOW POSITION/SIZE CHANGE
    def _reshape_window(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        box_diag = pow(pow(self.system.box_l[0], 2) + pow(
            self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)
        gluPerspective(
            40, 1.0 * w / h, self.specs['close_cut_distance'], self.specs['far_cut_distance'] * box_diag)
        glMatrixMode(GL_MODELVIEW)
        self.specs['window_size'][0] = w
        self.specs['window_size'][1] = h

    # INITS FOR GLUT FUNCTIONS
    def _init_openGL_callbacks(self):
        # OpenGl Callbacks
        def display():
            if self.hasParticleData and self.glutMainLoop_started:
                self._display_all()
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

        def redraw_on_idle():
            # DONT REPOST FASTER THAN 60 FPS
            self.draw_elapsed += (time.time() - self.draw_timer)
            if self.draw_elapsed > 1.0 / 60.0:
                self.draw_elapsed = 0
                glutPostRedisplay()
            self.draw_timer = time.time()
            return

        def reshape_callback(w, h):
            self._reshape_window(w, h)

        def close_window():
            os._exit(1)

        glutDisplayFunc(display)
        glutMouseFunc(mouse)
        glutKeyboardFunc(keyboard_down)
        glutKeyboardUpFunc(keyboard_up)
        glutSpecialFunc(keyboard_down)
        glutSpecialUpFunc(keyboard_up)
        glutReshapeFunc(reshape_callback)
        glutMotionFunc(motion)
        glutWMCloseFunc(close_window)

        glutIdleFunc(redraw_on_idle)

    def _init_timers(self):

        # TIMERS FOR register_callback
        def dummy_timer(index):
            self.timers[index][1]()
            glutTimerFunc(self.timers[index][0], dummy_timer, index)

        index = 0
        for t in self.timers:
            glutTimerFunc(t[0], dummy_timer, index)
            index += 1

        # HANDLE INPUT WITH 60FPS
        def timed_handle_input(data):
            self.keyboardManager.handle_input()
            glutTimerFunc(17, timed_handle_input, -1)

        glutTimerFunc(17, timed_handle_input, -1)

    # CLICKED ON PARTICLE: DRAG; CLICKED ON BACKGROUND: CAMERA
    def _mouse_motion(self, mousePos, mousePosOld, mouseButtonState):

        if self.specs['drag_enabled'] and self.dragId != -1:
            ppos = self.particles['pos'][self.dragId]
            viewport = glGetIntegerv(GL_VIEWPORT)
            mouseWorld = gluUnProject(
                mousePos[0], viewport[3] - mousePos[1], self.depth)

            self.dragExtForce = self.specs['drag_force'] * \
                (np.asarray(mouseWorld) - np.array(ppos))
            self.triggerSetParticleDrag = True
        else:
            self.camera.rotate_camera(mousePos, mousePosOld, mouseButtonState)

    # DRAW SCENE AGAIN WITHOUT LIGHT TO IDENTIFY PARTICLE ID BY PIXEL COLOR
    def _get_particle_id(self, pos, pos_old):

        glClearColor(0.0, 0.0, 0.0, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glLoadMatrixf(self.camera.modelview)

        glDisable(GL_LIGHTING)
        glDisable(GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            glDisable(GL_LIGHT1)
        self._draw_system_particles(colorById=True)
        viewport = glGetIntegerv(GL_VIEWPORT)

        readPixel = glReadPixelsui(
            pos[0], viewport[3] - pos[1], 1, 1, GL_RGB, GL_FLOAT)[0][0]
        depth = glReadPixelsf(
            pos[0], viewport[3] - pos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]

        pid = self._fcolor_to_id(readPixel)

        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            glEnable(GL_LIGHT1)
        glClearColor(self.specs['background_color'][0],
                     self.specs['background_color'][1],
                     self.specs['background_color'][2], 1.)

        return pid, depth

    def _id_to_fcolor(self, pid):
        pid += 1
        return [int(pid / (256 * 256)) / 255.0, int((pid % (256 * 256)) / 256) / 255.0, (pid % 256) / 255.0, 1.0]

    def _fcolor_to_id(self, fcol):
        if (fcol == [0, 0, 0]).all():
            return -1
        else:
            return 256 * 256 * int(fcol[0] * 255) + 256 * int(fcol[1] * 255) + int(fcol[2] * 255) - 1

    def _set_particle_drag(self, pos, pos_old):
        pid, depth = self._get_particle_id(pos, pos_old)
        self.dragId = pid

        if pid != -1:
            self.dragPosInitial = self.particles['pos'][self.dragId]
            self.extForceOld = self.particles['ext_force'][self.dragId][:]
            self.depth = depth

    def _reset_particle_drag(self, pos, pos_old):
        if self.dragId != -1:
            self.triggerResetParticleDrag = True

    def _get_particle_info(self, pos, pos_old):
        pid, depth = self._get_particle_id(pos, pos_old)
        if self.show_system_info:
            self.show_system_info = False
        elif pid == -1 and self.infoId == -1:
            self.show_system_info = True
            self.update_system_info()
        self.infoId = pid

    def _next_particle_info(self):
        self.infoId = (self.infoId + 1) % len(self.particles['pos'])

    def _previous_particle_info(self):
        self.infoId = (self.infoId - 1) % len(self.particles['pos'])

    # ESPRESSO RELATED INITS
    def _init_espresso_visualization(self):
        self.maxq = 0
        self.minq = 0

        self.dragId = -1
        self.infoId = -1
        self.show_system_info = False
        self.dragPosInitial = []
        self.extForceOld = []
        self.dragExtForceOld = []
        self.triggerResetParticleDrag = False
        self.triggerSetParticleDrag = False

        self.depth = 0

        self.imPos = [np.array([self.system.box_l[0], 0, 0]), np.array(
            [0, self.system.box_l[1], 0]), np.array([0, 0, self.system.box_l[2]])]

        # LOOK FOR LB ACTOR
        if self.specs['LB_draw_velocity_plane'] or self.specs['LB_draw_nodes'] or self.specs['LB_draw_node_boundaries']:
            for a in self.system.actors:
                types = []
                IF LB:
                    types.append(espressomd.lb.LBFluid)
                IF LB_GPU:
                    types.append(espressomd.lb.LBFluidGPU)

                # if type(a) == espressomd.lb.LBFluidGPU or type(a) == espressomd.lb.LBFluid
                if type(a) in types:
                    # if 'agrid' in pa:
                    self.lb_params = a.get_params()
                    self.lb = a
                    self.lb_is_cpu = type(a) == espressomd.lb.LBFluid
                    break

        if self.specs['LB_draw_velocity_plane']:

            if self.specs['LB_plane_axis'] == 0:
                pn = [1.0, 0.0, 0.0]
                self.lb_plane_b1 = [0.0, 1.0, 0.0]
                self.lb_plane_b2 = [0.0, 0.0, 1.0]
            elif self.specs['LB_plane_axis'] == 1:
                pn = [0.0, 1.0, 0.0]
                self.lb_plane_b1 = [1.0, 0.0, 0.0]
                self.lb_plane_b2 = [0.0, 0.0, 1.0]
            else:
                pn = [0.0, 0.0, 1.0]
                self.lb_plane_b1 = [1.0, 0.0, 0.0]
                self.lb_plane_b2 = [0.0, 1.0, 0.0]

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

            self._update_lb_velocity_plane()

        self._box_size_dependence()

    # BOX PLANES (NORMAL, ORIGIN) FOR PERIODIC BONDS
    def _box_size_dependence(self):

        if self.specs['draw_cells'] or self.specs['draw_nodes']:
            self._update_nodes()
            if self.specs['draw_cells']:
                self._update_cells()

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

        self.camera.center = np.array(self.system.box_l) * 0.5
        self.camera.update_modelview()

    # DEFAULT CONTROLS
    def _init_controls(self):

        # MOUSE LOOK/ROTATE/DRAG
        self.mouseManager.register_button(MouseButtonEvent(
            None, MouseFireEvent.FreeMotion, self._mouse_motion))

        self.mouseManager.register_button(MouseButtonEvent(
            3, MouseFireEvent.ButtonPressed, self.camera.move_backward))

        self.mouseManager.register_button(MouseButtonEvent(
            4, MouseFireEvent.ButtonPressed, self.camera.move_forward))

        # START/STOP DRAG
        if self.specs['drag_enabled']:
            self.mouseManager.register_button(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonPressed, self._set_particle_drag, True))
            self.mouseManager.register_button(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonReleased, self._reset_particle_drag, True))

        # PARTICLE INFORMATION
        self.mouseManager.register_button(MouseButtonEvent(
            GLUT_LEFT_BUTTON, MouseFireEvent.DoubleClick, self._get_particle_info, True))

        # CYCLE THROUGH PARTICLES
        self.keyboardManager.register_button(KeyboardButtonEvent(
            GLUT_KEY_LEFT, KeyboardFireEvent.Pressed, self._previous_particle_info, False))
        self.keyboardManager.register_button(KeyboardButtonEvent(
            GLUT_KEY_RIGHT, KeyboardFireEvent.Pressed, self._next_particle_info, False))

        # <SPACE> PAUSE INTEGRATION THREAD
        self.keyboardManager.register_button(KeyboardButtonEvent(
            ' ', KeyboardFireEvent.Pressed, self._pause, True))

        # <RETURN> TAKE SCREENSHOT
        self.keyboardManager.register_button(KeyboardButtonEvent(
            '\x0d', KeyboardFireEvent.Pressed, self._trigger_screenshot, True))

        # <ESCAPE> QUIT
        self.keyboardManager.register_button(KeyboardButtonEvent(
            '\x1b', KeyboardFireEvent.Pressed, self._quit, True))

        # CAMERA CONTROL VIA KEYBOARD
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
    def _quit(self):
        self.quit_savely = True

    def _pause(self):
        self.paused = not self.paused

    def _trigger_screenshot(self):
        self.take_screenshot = True

    # ASYNCHRONOUS PARALLEL CALLS OF glLight CAUSES SEG FAULTS, SO ONLY CHANGE
    # LIGHT AT CENTRAL display METHOD AND TRIGGER CHANGES
    def _set_camera_spotlight(self):
        if self.specs['spotlight_enabled']:
            p = self.camera.camPos
            fp = [p[0], p[1], p[2], 1]
            glLightfv(GL_LIGHT1, GL_POSITION, fp)
            glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, self.camera.state_target)

    def _init_camera(self):
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

        self.camera.set_camera(camPos=np.array(cp), camTarget=ct, camRight=cr, moveSpeed=0.5 *
                               box_diag / 17.0,  center=box_center)
        self._set_camera_spotlight()

    def _init_opengl(self):
        glutInit(self.specs['name'])
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)

        glutInitWindowSize(self.specs['window_size'][
                           0], self.specs['window_size'][1])

        glutCreateWindow(b"ESPResSo visualization")

        glClearColor(self.specs['background_color'][0], self.specs[
                     'background_color'][1], self.specs['background_color'][2], 1.)

        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_LINE_SMOOTH)
        glHint(GL_LINE_SMOOTH_HINT, GL_NICEST)

        # BAD FOR TRANSPARENT PARTICLES
        # glEnable(GL_CULL_FACE)
        # glCullFace(GL_BACK)

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
            glLightfv(GL_LIGHT0, GL_POSITION, np.array(
                self.specs['light_pos']).tolist())
        else:
            glLightfv(GL_LIGHT0, GL_POSITION,
                      (np.array(self.system.box_l) * 1.1).tolist())

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

        self.dl_sphere = glGenLists(1)


# END OF MAIN CLASS

# OPENGL DRAW WRAPPERS

def set_solid_material(color, material=[0.6, 1.0, 0.1, 0.4, 1.0]):
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  [
                 color[0] * material[0], color[1] * material[0], color[2] * material[0], material[4]])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  [
                 color[0] * material[1], color[1] * material[1], color[2] * material[1], material[4]])
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [
                 color[0] * material[2], color[1] * material[2], color[2] * material[2], material[4]])
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, int(material[3] * 128))


def draw_box(p0, s, color, material, width):
    glLineWidth(width)
    set_solid_material(color, material)
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

    glPopMatrix()


def draw_sphere(pos, radius, color, material, quality):
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    set_solid_material(color, material)
    glutSolidSphere(radius, quality, quality)
    glPopMatrix()


def draw_plane(edges, color, material):

    set_solid_material(color, material)

    glBegin(GL_QUADS)
    for e in edges:
        glVertex3f(e[0], e[1], e[2])
    glEnd()


def draw_points(points, pointsize, color, material):
    set_solid_material(color, material)
    glPointSize(pointsize)
    glBegin(GL_POINTS)
    for p in points:
        glVertex3f(p[0], p[1], p[2])
    glEnd()


def draw_cylinder(posA, posB, radius, color, material, quality, draw_caps=False):
    set_solid_material(color, material)
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
    # the rotation axis is the cross product between z and d
    vz = np.cross([0.0, 0.0, 1.0], d)
    # get the angle using a dot product
    angle = 180.0 / np.pi * acos(d[2] / np.linalg.norm(d))

    return angle, vz[0], vz[1]


def draw_ellipsoid(pos, semiaxis_a, semiaxis_b, semiaxis_c, color, material, quality):
    set_solid_material(color, material)
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glScalef(semiaxis_a, semiaxis_b, semiaxis_c)
    glutSolidSphere(1, quality, quality)
    glPopMatrix()


def get_extra_clip_plane():

    # ON SOME HARDWARE (e.g. MAC) only 6 CLIP PLANES ARE ALLOWED,
    # SO CAPPING OF BOX BOUNDARIES AND ADDITIONAL SHAPE CLIP PLANES
    # ARE NOT POSSIBLE. THIS WILL CAUSE THE SHAPES THAT NEED ADDITIONAL
    # CLIP PLANES TO NOT BE CLIPPED ON ONE FACE OF THE BOX

    if GL_MAX_CLIP_PLANES > 6:
        return GL_CLIP_PLANE0 + 6
    else:
        return GL_CLIP_PLANE0


def draw_simple_pore(center, axis, length, radius, smoothing_radius, max_box_l, color, material, quality):

    clip_plane = get_extra_clip_plane()

    set_solid_material(color, material)
    glPushMatrix()
    quadric = gluNewQuadric()

    # basic position and orientation
    glTranslate(center[0], center[1], center[2])
    ax, rx, ry = rotation_helper(axis)
    glRotatef(ax, rx, ry, 0.0)
    # cylinder
    glTranslate(0, 0, -0.5 * length + smoothing_radius)
    gluCylinder(quadric, radius, radius, length - 2 *
                smoothing_radius, quality, quality)
    # torus segment

    glEnable(clip_plane)
    glClipPlane(clip_plane, (0, 0, -1, 0))
    glutSolidTorus(smoothing_radius, (radius +
                                      smoothing_radius), quality, quality)
    glDisable(clip_plane)
    # wall
    glTranslate(0, 0, -smoothing_radius)
    gluPartialDisk(quadric, radius + smoothing_radius,
                   2.0 * max_box_l, quality, 1, 0, 360)
    # torus segment
    glTranslate(0, 0, length - smoothing_radius)
    glEnable(clip_plane)
    glClipPlane(clip_plane, (0, 0, 1, 0))
    glutSolidTorus(smoothing_radius, (radius +
                                      smoothing_radius), quality, quality)
    glDisable(clip_plane)
    # wall
    glTranslate(0, 0, smoothing_radius)
    gluPartialDisk(quadric, radius + smoothing_radius,
                   2.0 * max_box_l, quality, 1, 0, 360)

    glPopMatrix()


def draw_sphero_cylinder(posA, posB, radius, color, material, quality):
    set_solid_material(color, material)
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
    clip_plane = get_extra_clip_plane()
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
    DoubleClick = 4


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
        self.mouseEventsDoubleClick = []
        self.mouseState = {}
        self.mouseState[GLUT_LEFT_BUTTON] = GLUT_UP
        self.mouseState[GLUT_MIDDLE_BUTTON] = GLUT_UP
        self.mouseState[GLUT_RIGHT_BUTTON] = GLUT_UP
        self.mouseState['3'] = GLUT_UP  # WHEEL
        self.mouseState['4'] = GLUT_UP  # WHEEL
        self.pressedTime = {}
        self.pressedTime[GLUT_LEFT_BUTTON] = 0
        self.pressedTime[GLUT_MIDDLE_BUTTON] = 0
        self.pressedTime[GLUT_RIGHT_BUTTON] = 0
        self.pressedTime[3] = 0
        self.pressedTime[4] = 0
        self.pressedTimeOld = {}
        self.pressedTimeOld[GLUT_LEFT_BUTTON] = 0
        self.pressedTimeOld[GLUT_MIDDLE_BUTTON] = 0
        self.pressedTimeOld[GLUT_RIGHT_BUTTON] = 0
        self.pressedTimeOld[3] = 0
        self.pressedTimeOld[4] = 0

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
        elif mouseEvent.fireEvent == MouseFireEvent.DoubleClick:
            self.mouseEventsDoubleClick.append(mouseEvent)

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

        if state == GLUT_DOWN:
            self.pressedTimeOld[button] = self.pressedTime[button]
            self.pressedTime[button] = time.time()

        for me in self.mouseEventsDoubleClick:
            if me.button == button and state == GLUT_DOWN and self.pressedTime[button] - self.pressedTimeOld[button] < 0.25:
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

    def __init__(self):
        pass

    def set_camera(self, camPos=np.array([0, 0, 1]), camTarget=np.array([0, 0, 0]), camRight=np.array([1.0, 0.0, 0.0]), moveSpeed=0.5, rotSpeed=0.001, globalRotSpeed=3.0, center=np.array([0, 0, 0])):
        self.moveSpeed = moveSpeed
        self.lookSpeed = rotSpeed
        self.globalRotSpeed = globalRotSpeed

        self.center = center

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
        self.rotate_system_y(0.01 * self.globalRotSpeed)

    def rotate_system_XR(self):
        self.rotate_system_y(-0.01 * self.globalRotSpeed)

    def rotate_system_YL(self):
        self.rotate_system_z(0.01 * self.globalRotSpeed)

    def rotate_system_YR(self):
        self.rotate_system_z(-0.01 * self.globalRotSpeed)

    def rotate_system_ZL(self):
        self.rotate_system_x(0.01 * self.globalRotSpeed)

    def rotate_system_ZR(self):
        self.rotate_system_x(-0.01 * self.globalRotSpeed)

    def rotate_camera(self, mousePos, mousePosOld, mouseButtonState):
        dm = mousePos - mousePosOld

        if mouseButtonState[GLUT_LEFT_BUTTON] == GLUT_DOWN:
            if dm[0] != 0:
                self.rotate_system_y(dm[0] * 0.001 * self.globalRotSpeed)
            if dm[1] != 0:
                self.rotate_system_x(dm[1] * 0.001 * self.globalRotSpeed)
        elif mouseButtonState[GLUT_RIGHT_BUTTON] == GLUT_DOWN:
            self.state_pos[0] -= 0.05 * dm[0] * self.moveSpeed
            self.state_pos[1] += 0.05 * dm[1] * self.moveSpeed
            self.update_modelview()
        elif mouseButtonState[GLUT_MIDDLE_BUTTON] == GLUT_DOWN:
            self.state_pos[2] += 0.05 * dm[1] * self.moveSpeed
            self.rotate_system_z(dm[0] * 0.001 * self.globalRotSpeed)

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

    def rotate_system_z(self, da):
        self.rotate_vector(self.state_right, da, self.state_target)
        self.rotate_vector(self.state_up, da, self.state_target)
        self.update_modelview()

    def rotate_system_x(self, da):
        self.rotate_vector(self.state_target, da, self.state_right)
        self.state_up = np.cross(self.state_right, self.state_target)
        self.update_modelview()

    def rotate_system_y(self, da):
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
        x = (self[0] * q[3]) + (self[3] * q[0]) + \
            (self[1] * q[2]) - (self[2] * q[1])
        y = (self[1] * q[3]) + (self[3] * q[1]) + (
            self[2] * q[0]) - (self[0] * q[2])
        z = (self[2] * q[3]) + (self[3] * q[2]) + (
            self[0] * q[1]) - (self[1] * q[0])
        return Quaternion(x, y, z, w)
