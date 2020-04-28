# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import ctypes
import math
import os
import sys
import time
from threading import Thread

import OpenGL.GL
import OpenGL.GLE
import OpenGL.GLU
import OpenGL.GLUT
import numpy as np
from matplotlib.pyplot import imsave

import espressomd
from .particle_data import ParticleHandle


class openGLLive():
    """
    This class provides live visualization using pyOpenGL.
    Use the update method to push your current simulation state after
    integrating. Modify the appearance with a list of keywords.
    Timed callbacks can be registered via the :meth:`register_callback` method
    and keyboard callbacks via :meth:`keyboard_manager.register_button()
    <espressomd.visualization_opengl.keyboard_manager.register_button>`.

    Parameters
    ----------

    system : :class:`espressomd.system.System`
    window_size : (2,) array_like of :obj:`int`, optional
        Size of the visualizer window in pixels.
    name : :obj:`str`, optional
        The name of the visualizer window.
    background_color : (3,) array_like of :obj:`int`, optional
        RGB of the background.
    periodic_images : (3,) array_like of :obj:`int`, optional
        Periodic repetitions on both sides of the box in xyz-direction.
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
    camera_position : :obj:`str` or (3,) array_like of :obj:`float`, optional
        Initial camera position. Use ``'auto'`` (default) for shifted position in z-direction.
    camera_target : :obj:`str` or (3,) array_like of :obj:`float`, optional
        Initial camera target. Use ``'auto'`` (default) to look towards the system center.
    camera_right : (3,) array_like of :obj:`float`, optional
        Camera right vector in system coordinates. Default is ``[1, 0, 0]``
    particle_sizes : :obj:`str` or array_like :obj:`float` or callable, optional
        * ``'auto'`` (default): The Lennard-Jones sigma value of the
          self-interaction is used for the particle diameter.
        * callable: A lambda function with one argument. Internally,
          the numerical particle type is passed to the lambda
          function to determine the particle radius.
        * list: A list of particle radii, indexed by the particle type.
    particle_coloring : :obj:`str`, optional
        * ``'auto'`` (default): Colors of charged particles are
          specified by ``particle_charge_colors``, neutral particles
          by ``particle_type_colors``.
        * ``'charge'``: Minimum and maximum charge of all particles is determined by the
          visualizer. All particles are colored by a linear
          interpolation of the two colors given by
          ``particle_charge_colors`` according to their charge.
        * ``'type'``: Particle colors are specified by ``particle_type_colors``,
          indexed by their numerical particle type.
        * ``'node'``: Color according to the node the particle is on.
    particle_type_colors : array_like :obj:`float`, optional
        Colors for particle types.
    particle_type_materials : array_like :obj:`str`, optional
        Materials of the particle types.
    particle_charge_colors : (2,) array_like of :obj:`float`, optional
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
        Colors of the constraints by type.
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
        Materials of ext_force arrows for different particle types.
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
    light_pos : (3,) array_like of :obj:`float`, optional
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

    Notes
    -----
    The visualization of some constraints is either improved by or even relies
    on the presence of an installed OpenGL Extrusion library on your system.

    """

    # MATERIALS
    materials = {
        'bright': [0.9, 1.0, 0.8, 0.4, 1.0],
        'medium': [0.6, 0.8, 0.2, 0.4, 1.0],
        'dark': [0.4, 0.5, 0.1, 0.4, 1.0],
        'transparent1': [0.6, 0.8, 0.2, 0.5, 0.8],
        'transparent2': [0.6, 0.8, 0.2, 0.5, 0.4],
        'transparent3': [0.6, 0.8, 0.2, 0.5, 0.2],
        'rubber': [0, 0.4, 0.7, 0.078125, 1.0],
        'chrome': [0.25, 0.4, 0.774597, 0.6, 1.0],
        'plastic': [0, 0.55, 0.7, 0.25, 1.0],
        'steel': [0.25, 0.38, 0, 0.32, 1.0]
    }

    def __init__(self, system, **kwargs):
        # DEFAULT PROPERTIES
        self.specs = {
            'window_size': [800, 800],
            'name': 'ESPResSo Visualization',

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
            'particle_type_colors':
                [[1, 1, 0], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 1], [1, 0.5, 0], [0.5, 0, 1]],
            'particle_type_materials': ['medium'],
            'particle_charge_colors': [[1, 0, 0], [0, 1, 0]],

            'draw_constraints': True,
            'rasterize_pointsize': 10,
            'rasterize_resolution': 75.0,
            'constraint_type_colors':
                [[0.5, 0.5, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5],
                    [0.5, 0.5, 0], [0, 0, 0.5], [0.5, 0, 0]],
            'constraint_type_materials': ['transparent1'],

            'draw_bonds': True,
            'bond_type_radius': [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            'bond_type_colors':
                [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'bond_type_materials': ['medium'],

            'ext_force_arrows': False,
            'ext_force_arrows_type_scale': [1.0],
            'ext_force_arrows_type_colors':
                [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'ext_force_arrows_type_materials': ['transparent2'],
            'ext_force_arrows_type_radii': [0.2],

            'velocity_arrows': False,
            'velocity_arrows_type_scale': [1.0],
            'velocity_arrows_type_colors':
                [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'velocity_arrows_type_materials': ['transparent2'],
            'velocity_arrows_type_radii': [0.2],

            'force_arrows': False,
            'force_arrows_type_scale': [1.0],
            'force_arrows_type_colors':
                [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
            'force_arrows_type_materials': ['transparent2'],
            'force_arrows_type_radii': [0.2],

            'director_arrows': False,
            'director_arrows_type_scale': [1.0],
            'director_arrows_type_colors':
                [[1, 1, 1], [1, 0, 1], [0, 0, 1], [0, 1, 1],
                    [1, 1, 0], [1, 0.5, 0], [0.5, 0, 1]],
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
            'light_colors':
                [[0.1, 0.1, 0.2], [0.9, 0.9, 0.9], [1.0, 1.0, 1.0]],
            'light_brightness': 1.0,
            'light_size': 'auto',

            'spotlight_enabled': False,
            'spotlight_colors':
                [[0.2, 0.2, 0.3], [0.5, 0.5, 0.5], [1.0, 1.0, 1.0]],
            'spotlight_angle': 45,
            'spotlight_focus': 1,
            'spotlight_brightness': 0.6,

            'drag_enabled': False,
            'drag_force': 3.0
        }

        # OVERWRITE WITH USER PROPERTIES
        for key in kwargs:
            if key not in self.specs:
                raise ValueError(f'{key} is no valid visualization property')
            else:
                self.specs[key] = kwargs[key]

        # DEPENDENCIES
        if not espressomd.has_features('EXTERNAL_FORCES'):
            self.specs['drag_enabled'] = False
            self.specs['ext_force_arrows'] = False

        if not espressomd.has_features('ROTATION'):
            self.specs['director_arrows'] = False

        if not espressomd.has_features('LB_BOUNDARIES') and \
                not espressomd.has_features('LB_BOUNDARIES_GPU'):
            self.specs['LB_draw_boundaries'] = False
            self.specs['LB_draw_node_boundaries'] = False

        # ESPRESSO RELATED INITS THAT ARE KNOWN ONLY WHEN RUNNING THE
        # INTEGRATION LOOP ARE CALLED ONCE IN UPDATE LOOP:
        # CONSTRAINTS, NODE BOXES, CELL BOXES, CHARGE RANGE, BONDS

        # ESPRESSO RELATED INITS THAT ARE KNOWN ON INSTANTIATION GO HERE:

        # CONTENT OF PARTICLE DATA
        self.has_particle_data = {
            'velocity': self.specs['velocity_arrows'],
            'force': self.specs['force_arrows'],
            'ext_force': self.specs['ext_force_arrows'] or
            self.specs['drag_enabled']}

        if espressomd.has_features('ELECTROSTATICS'):
            self.has_particle_data['charge'] = \
                self.specs['particle_coloring'] == 'auto' or \
                self.specs['particle_coloring'] == 'charge'
        else:
            self.has_particle_data['charge'] = False
        self.has_particle_data['director'] = self.specs['director_arrows']
        self.has_particle_data['node'] = self.specs['particle_coloring'] == 'node'

        # PARTICLE INFO OF HIGHLIGHTED PARTICLE: COLLECT PARTICLE ATTRIBUTES
        self.highlighted_particle = {}
        self.particle_attributes = []
        for d in dir(ParticleHandle):
            if isinstance(getattr(ParticleHandle, d),
                          type(ParticleHandle.pos)):
                if d not in ["pos_folded"]:
                    self.particle_attributes.append(d)
        self.max_len_attr = max([len(a) for a in self.particle_attributes])

        # FIXED COLORS FROM INVERSE BACKGROUND COLOR FOR GOOD CONTRAST
        self.inverse_bg_color = \
            np.array([1 - self.specs['background_color'][0],
                      1 - self.specs['background_color'][1],
                      1 - self.specs['background_color'][2], 1.0])

        self.node_box_color = np.copy(self.inverse_bg_color)
        self.node_box_color[0] += 0.5 * (0.5 - self.node_box_color[0])

        self.cell_box_color = np.copy(self.inverse_bg_color)
        self.cell_box_color[1] += 0.5 * (0.5 - self.cell_box_color[1])

        self.lb_box_color = np.copy(self.inverse_bg_color)
        self.lb_box_color[2] = 0.5

        self.lb_box_color_boundary = np.copy(self.inverse_bg_color)
        self.lb_box_color_boundary[1] = 0.5

        self.text_color = np.copy(self.inverse_bg_color)

        # INCREASE LINE THICKNESS IF NODE/CELL BOX IS ENABLED
        self.line_width_fac = 1.0
        if self.specs['draw_nodes']:
            self.line_width_fac += 0.5
        if self.specs['draw_cells']:
            self.line_width_fac += 0.5

        # HAS PERIODIC IMAGES
        self.has_images = any(i != 0 for i in self.specs['periodic_images'])
        self.image_vectors = [
            list(range(-i, i + 1)) for i in self.specs['periodic_images']]

        # INITS
        self.system = system
        self.system_info = {}

        self.last_T = -1
        self.last_box_l = self.system.box_l
        self.fps_last = 0
        self.fps = 0
        self.fps_count = 0

        self.glut_main_loop_started = False
        self.screenshot_initialized = False
        self.hasParticleData = False
        self.quit_safely = False
        self.paused = False
        self.take_screenshot = False
        self.screenshot_captured = False

        self.keyboard_manager = KeyboardManager()
        self.mouse_manager = MouseManager()
        self.camera = self._init_camera()

        self.timers = []
        self.particles = {}

        self.update_elapsed = 0
        self.update_timer = 0
        self.draw_elapsed = 0
        self.draw_timer = 0

        self.trigger_set_particle_drag = False
        self.trigger_reset_particle_drag = False
        self.drag_id = -1

        # LIST OF [[px,py],[string]] FOR USER DEFINED TEXT
        self.user_texts = []

    def update_system_info(self):
        """Update the information stored in dict self.system_info.

        The dictionary is used for showing system information, such as particle
        or constraint information, in the visualizer window.

        """
        # SYSTEM INFORMATION
        self.system_info = {
            'Actors': [], 'Non-bonded interactions': [],
            'Bonded interactions': [b for b in self.system.bonded_inter],
            'Constraints': [], 'Thermostat': self.system.thermostat.get_state()
        }

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

        # COLLECT ALL ACTIVE NON-BONDED INTERACTIONS
        all_non_bonded_inters = [
            x for x in dir(self.system.non_bonded_inter[0, 0])
            if not x.startswith('__') and x != 'type1' and x != 'type2']
        for t1 in all_types:
            for t2 in all_types:
                for check_nb in all_non_bonded_inters:
                    nb = getattr(
                        self.system.non_bonded_inter[t1, t2], check_nb)
                    if nb is not None and nb.is_active():
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
        """Renders the current state into an image file at ``path`` with
        dimensions of ``specs['window_size']`` in PNG format.

        """
        # ON FIRST CALL: INIT AND CREATE BUFFERS
        if not self.screenshot_initialized:
            self.screenshot_initialized = True
            self._init_opengl()

            # CREATE BUFFERS THAT CAN BE LARGER THAN THE SCREEN FRAME BUFFER
            fbo = OpenGL.GL.glGenFramebuffers(1)
            OpenGL.GL.glBindFramebuffer(OpenGL.GL.GL_FRAMEBUFFER, fbo)
            # COLOR BUFFER
            rbo = OpenGL.GL.glGenRenderbuffers(1)
            OpenGL.GL.glBindRenderbuffer(OpenGL.GL.GL_RENDERBUFFER, rbo)
            OpenGL.GL.glRenderbufferStorage(
                OpenGL.GL.GL_RENDERBUFFER,
                OpenGL.GL.GL_RGB,
                self.specs['window_size'][0],
                self.specs['window_size'][1])
            OpenGL.GL.glFramebufferRenderbuffer(
                OpenGL.GL.GL_FRAMEBUFFER,
                OpenGL.GL.GL_COLOR_ATTACHMENT0,
                OpenGL.GL.GL_RENDERBUFFER,
                rbo)
            # DEPTH BUFFER
            dbo = OpenGL.GL.glGenRenderbuffers(1)
            OpenGL.GL.glBindRenderbuffer(OpenGL.GL.GL_RENDERBUFFER, dbo)
            OpenGL.GL.glRenderbufferStorage(
                OpenGL.GL.GL_RENDERBUFFER, OpenGL.GL.GL_DEPTH_COMPONENT,
                self.specs['window_size'][0], self.specs['window_size'][1])
            OpenGL.GL.glFramebufferRenderbuffer(
                OpenGL.GL.GL_FRAMEBUFFER,
                OpenGL.GL.GL_DEPTH_ATTACHMENT,
                OpenGL.GL.GL_RENDERBUFFER,
                dbo)

            self._reshape_window(
                self.specs['window_size'][0], self.specs['window_size'][1])
            OpenGL.GLUT.glutHideWindow()

        # INIT AND UPDATE ESPRESSO
        self._init_espresso_visualization()
        self._initial_espresso_updates()

        # DRAW
        OpenGL.GL.glClear(
            OpenGL.GL.GL_COLOR_BUFFER_BIT | OpenGL.GL.GL_DEPTH_BUFFER_BIT)
        OpenGL.GL.glLoadMatrixf(self.camera.modelview)
        self._draw_system()

        # READ THE PIXELS
        OpenGL.GL.glReadBuffer(OpenGL.GL.GL_COLOR_ATTACHMENT0)
        data = OpenGL.GL.glReadPixels(
            0,
            0,
            self.specs['window_size'][0],
            self.specs['window_size'][1],
            OpenGL.GL.GL_RGB,
            OpenGL.GL.GL_FLOAT)

        # RESHAPE THE DATA
        data = np.flipud(data.reshape((data.shape[1], data.shape[0], 3)))

        # SAVE TO IMAGE
        imsave(path, data)

    def run(self, integration_steps=1):
        """Convenience method with a simple integration thread.

        """

        def main():
            while True:

                self.update()

                if self.paused:
                    time.sleep(0.0001)  # sleep(0) is worse
                else:
                    try:
                        self.system.integrator.run(integration_steps)
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
        self._init_OpenGL_callbacks()
        self._init_timers()

        # START THE BLOCKING MAIN LOOP
        self.glut_main_loop_started = True
        OpenGL.GLUT.glutMainLoop()

    def update(self):
        """Update method to be called after integration.
        Changes of ESPResSo system can only happen here.

        """
        if self.glut_main_loop_started:

            # UPDATE ON STARTUP
            if not self.hasParticleData:
                self._initial_espresso_updates()
                self.hasParticleData = True

            # UPDATES
            self.update_elapsed += (time.time() - self.update_timer)
            if self.update_elapsed > 1.0 / self.specs['update_fps']:
                self.update_elapsed = 0

                # ES UPDATES WHEN SYSTEM HAS PROPAGATED. ALSO UPDATE ON PAUSE
                # FOR PARTICLE INFO
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
                for c in self.keyboard_manager.userCallbackStack:
                    c()
                self.keyboard_manager.userCallbackStack = []

            self.update_timer = time.time()

            # DRAG PARTICLES
            if self.specs['drag_enabled']:
                if self.trigger_set_particle_drag and self.drag_id != -1:
                    self.system.part[self.drag_id].ext_force = self.dragExtForce
                    self.trigger_set_particle_drag = False
                elif self.trigger_reset_particle_drag and self.drag_id != -1:
                    self.system.part[self.drag_id].ext_force = self.extForceOld
                    self.trigger_reset_particle_drag = False
                    self.drag_id = -1

        # Escape was pressed: wait for ES to finish, then call sys exit from
        # main thread
        if self.quit_safely:
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

        if self.info_id != -1:
            for attr in self.particle_attributes:
                self.highlighted_particle[attr] = getattr(
                    self.system.part[self.info_id], attr)

    def _update_lb_velocity_plane(self):
        if self.lb_is_cpu:
            self._update_lb_velocity_plane_cpu()
        else:
            self._update_lb_velocity_plane_gpu()

    def _update_lb_velocity_plane_cpu(self):
        agrid = self.lb_params['agrid']
        self.lb_plane_vel = []
        ng = self.specs['LB_plane_ngrid']
        for xi in range(ng):
            for xj in range(ng):
                pp = np.copy((self.lb_plane_p + xi * 1.0 / ng * self.lb_plane_b1 +
                              xj * 1.0 / ng * self.lb_plane_b2) % self.system.box_l)
                i, j, k = (int(ppp / agrid) for ppp in pp)
                lb_vel = np.copy(self.lb[i, j, k].velocity)
                self.lb_plane_vel.append([pp, lb_vel])

    def _update_lb_velocity_plane_gpu(self):
        ng = self.specs['LB_plane_ngrid']
        col_pos = []
        for xi in range(ng):
            for xj in range(ng):
                p = np.array((self.lb_plane_p + xi * 1.0 / ng * self.lb_plane_b1 +
                              xj * 1.0 / ng * self.lb_plane_b2) % self.system.box_l)
                col_pos.append(p)

        lb_vels = self.lb.get_interpolated_fluid_velocity_at_positions(
            np.array(col_pos))
        self.lb_plane_vel = []
        for p, v in zip(col_pos, lb_vels):
            self.lb_plane_vel.append([p, v])

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
        self.local_box_l = np.array([0, 0, 0])
        for i in range(self.system.cell_system.node_grid[0]):
            for j in range(self.system.cell_system.node_grid[1]):
                for k in range(self.system.cell_system.node_grid[2]):
                    self.node_box_origins.append(
                        np.array([i, j, k]) * self.local_box_l)

    # GET THE _update_constraints DATA
    def _update_constraints(self):
        """Collect shapes and interaction type (for coloring) from constraints
        and update the list self.shapes with instances of the respective Shape
        (or child) classes.

        """
        self.shapes = []

        # helper functions and shape dictionary
        def unpack_shapes(shapes_list):
            """Return a list of shapes, where all unions have been unpacked.

            """
            shapes = []
            try:
                # recursively unpack unions
                for shape in shapes_list:
                    shapes.extend(unpack_shapes(shape))
            except TypeError:
                # otherwise just append the shape
                shapes.append(shapes_list)
            return shapes

        def shape_arguments(shape, part_type):
            """Helper function returning argument list for initialization of
            Shape class (or respective child classes).

            """
            arguments = [
                shape, part_type,
                self._modulo_indexing(self.specs['constraint_type_colors'],
                                      part_type),
                self.materials[self._modulo_indexing(
                    self.specs['constraint_type_materials'], part_type)],
                self.specs['quality_constraints'], self.system.box_l,
                self.specs['rasterize_resolution'],
                self.specs['rasterize_pointsize']
            ]
            return arguments

        shape_mapping = {
            'Shapes::Cylinder': Cylinder,
            'Shapes::Ellipsoid': Ellipsoid,
            'Shapes::HollowConicalFrustum': HollowConicalFrustum,
            'Shapes::SimplePore': SimplePore,
            'Shapes::Slitpore': Slitpore,
            'Shapes::Sphere': Sphere,
            'Shapes::SpheroCylinder': Spherocylinder,
            'Shapes::Wall': Wall
        }

        # collect shapes using the helper functions
        for constraint in self.system.constraints:
            if isinstance(constraint,
                          espressomd.constraints.ShapeBasedConstraint):
                part_type = constraint.get_parameter('particle_type')
                shape = constraint.get_parameter('shape')
                for sub_shape in unpack_shapes(shape):
                    arguments = shape_arguments(sub_shape, part_type)
                    try:
                        self.shapes.append(
                            shape_mapping[sub_shape.name()](*arguments))
                    except KeyError:
                        self.shapes.append(Shape(*arguments))

        if self.specs['LB_draw_boundaries']:
            ni = 0
            for constraint in self.system.lbboundaries:
                if isinstance(constraint, espressomd.lbboundaries.LBBoundary):
                    part_type = ni
                    ni += 1
                    shape = constraint.get_parameter('shape')
                    for sub_shape in unpack_shapes(shape):
                        arguments = shape_arguments(sub_shape, part_type)
                        try:
                            self.shapes.append(
                                shape_mapping[sub_shape.name()](*arguments))
                        except KeyError:
                            self.shapes.append(Shape(*arguments))

    # GET THE BOND DATA, SO FAR CALLED ONCE UPON INITIALIZATION
    def _update_bonds(self):
        if self.specs['draw_bonds']:
            self.bonds = []
            for i, particle in enumerate(self.system.part):
                for bond in particle.bonds:
                    # b[0]: Bond, b[1:] Partners
                    bond_type = bond[0].type_number()
                    if len(bond) == 4:
                        self.bonds.append([i, bond[1], bond_type])
                        self.bonds.append([i, bond[2], bond_type])
                        self.bonds.append([bond[2], bond[3], bond_type])
                    else:
                        for bond_partner in bond[1:]:
                            self.bonds.append([i, bond_partner, bond_type])

    def _draw_text(self, x, y, text, color,
                   font=OpenGL.GLUT.GLUT_BITMAP_9_BY_15):
        OpenGL.GL.glColor(color)
        OpenGL.GL.glWindowPos2f(x, y)
        for ch in text:
            OpenGL.GLUT.glutBitmapCharacter(font, ctypes.c_int(ord(ch)))

    # DRAW CALLED AUTOMATICALLY FROM GLUT DISPLAY FUNC
    def _draw_system(self):
        if self.specs['LB_draw_velocity_plane']:
            self._draw_lb_vel()

        if self.specs['draw_axis']:
            axis_fac = 0.2
            axis_r = np.min(self.system.box_l) / 50.0
            draw_arrow([0, 0, 0], [self.system.box_l[0] * axis_fac, 0, 0], axis_r,
                       [1, 0, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, self.system.box_l[2] * axis_fac, 0], axis_r,
                       [0, 1, 0, 1], self.materials['chrome'], self.specs['quality_arrows'])
            draw_arrow([0, 0, 0], [0, 0, self.system.box_l[2] * axis_fac], axis_r,
                       [0, 0, 1, 1], self.materials['chrome'], self.specs['quality_arrows'])

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
        draw_box([0, 0, 0], self.system.box_l, self.inverse_bg_color,
                 self.materials['medium'], 2.0 * self.line_width_fac)

    def _draw_nodes(self):
        for n in self.node_box_origins:
            draw_box(n, self.local_box_l, self.node_box_color,
                     self.materials['transparent1'], 1.5 * self.line_width_fac)

    def _draw_cells(self):
        for n in self.node_box_origins:
            for c in self.cell_box_origins:
                draw_box(
                    c + n, self.cell_size, self.cell_box_color,
                    self.materials['transparent1'], 0.75 * self.line_width_fac)

    def _draw_lb_grid(self):
        a = self.lb_params['agrid']
        cell_size = np.array([a] * 3)
        dims = np.rint(np.array(self.system.box_l) / a)
        for i in range(int(dims[0])):
            for j in range(int(dims[1])):
                for k in range(int(dims[2])):
                    n = np.array([i, j, k]) * cell_size
                    if self.specs['LB_draw_node_boundaries'] \
                            and self.lb[i, j, k].boundary:
                        draw_box(n, cell_size, self.lb_box_color_boundary,
                                 self.materials['transparent2'], 5.0)
                    if self.specs['LB_draw_nodes'] \
                            and not self.lb[i, j, k].boundary:
                        draw_box(n, cell_size, self.lb_box_color,
                                 self.materials['transparent2'], 1.5)

    def _draw_constraints(self):

        # CLIP BORDERS OF SIMULATION BOX
        for i in range(6):
            OpenGL.GL.glEnable(OpenGL.GL.GL_CLIP_PLANE0 + i)
            OpenGL.GL.glClipPlane(OpenGL.GL.GL_CLIP_PLANE0 + i,
                                  self.box_eqn[i, :])

        for shape in self.shapes:
            shape.draw()

        for i in range(6):
            OpenGL.GL.glDisable(OpenGL.GL.GL_CLIP_PLANE0 + i)

    def _determine_radius(self, part_type):
        def radius_by_lj(part_type):
            try:
                radius = self.system.non_bonded_inter[part_type, part_type].lennard_jones.get_params()[
                    'sigma'] * 0.5
                if radius == 0:
                    radius = self.system.non_bonded_inter[part_type, part_type].wca.get_params()[
                        'sigma'] * 0.5
            except BaseException:
                radius = 0.5
            if radius == 0:
                radius = 0.5
            return radius

        if self.specs['particle_sizes'] == 'auto':
            radius = radius_by_lj(part_type)
        elif callable(self.specs['particle_sizes']):
            radius = self.specs['particle_sizes'](part_type)
        else:
            try:
                radius = self._modulo_indexing(
                    self.specs['particle_sizes'], part_type)
            except RuntimeError:
                radius = self.radius_by_lj(part_type)
        return radius

    def _draw_system_particles(self, color_by_id=False):
        part_ids = range(len(self.particles['pos']))
        part_type = -1
        reset_material = False

        for part_id in part_ids:
            part_type_last = part_type
            part_type = int(self.particles['type'][part_id])

            # Only change material if type/charge has changed, color_by_id or
            # material was reset by arrows
            if reset_material or color_by_id or not part_type == part_type_last or \
                    part_id == self.drag_id or part_id == self.info_id or self.specs['particle_coloring'] == 'node':
                reset_material = False

                radius = self._determine_radius(part_type)

                m = self._modulo_indexing(
                    self.specs['particle_type_materials'], part_type)
                material = self.materials[m]

                if color_by_id:
                    color = self._id_to_fcolor(part_id)
                    OpenGL.GL.glColor(color)
                else:
                    if self.specs['particle_coloring'] == 'auto':
                        # Color auto: Charge then Type
                        if self.has_particle_data['charge'] and self.particles['charge'][part_id] != 0:
                            color = self._color_by_charge(
                                self.particles['charge'][part_id])
                            reset_material = True
                        else:
                            color = self._modulo_indexing(
                                self.specs['particle_type_colors'], part_type)
                    elif self.specs['particle_coloring'] == 'charge':
                        color = self._color_by_charge(
                            self.particles['charge'][part_id])
                        reset_material = True
                    elif self.specs['particle_coloring'] == 'type':
                        color = self._modulo_indexing(
                            self.specs['particle_type_colors'], part_type)
                    elif self.specs['particle_coloring'] == 'node':
                        color = self._modulo_indexing(
                            self.specs['particle_type_colors'], self.particles['node'][part_id])

                    # Invert color of highlighted particle
                    if part_id == self.drag_id or part_id == self.info_id:
                        reset_material = True
                        color = [1 - color[0], 1 - color[1],
                                 1 - color[2]]

                    set_solid_material(color, material)

                # Create a new display list, used until next material/color
                # change
                OpenGL.GL.glNewList(self.dl_sphere, OpenGL.GL.GL_COMPILE)
                OpenGL.GLUT.glutSolidSphere(
                    radius, self.specs['quality_particles'], self.specs['quality_particles'])
                OpenGL.GL.glEndList()

            self._redraw_sphere(self.particles['pos'][part_id])

            if self.has_images:
                for imx in self.image_vectors[0]:
                    for imy in self.image_vectors[1]:
                        for imz in self.image_vectors[2]:
                            if imx != 0 or imy != 0 or imz != 0:
                                offset = [imx, imy, imz] * self.system.box_l
                                self._redraw_sphere(
                                    self.particles['pos'][part_id] + offset)

            if espressomd.has_features('EXTERNAL_FORCES'):
                if self.specs['ext_force_arrows'] or part_id == self.drag_id:
                    if any(
                            v != 0 for v in self.particles['ext_force'][part_id]):
                        if part_id == self.drag_id:
                            sc = 1
                        else:
                            sc = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_scale'], part_type)
                        if sc > 0:
                            arrow_col = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_colors'], part_type)
                            arrow_radius = self._modulo_indexing(
                                self.specs['ext_force_arrows_type_radii'], part_type)
                            draw_arrow(self.particles['pos'][part_id], np.array(
                                self.particles['ext_force'][part_id]) * sc, arrow_radius, arrow_col,
                                self.materials['chrome'], self.specs['quality_arrows'])
                            reset_material = True

            if self.specs['velocity_arrows']:
                self._draw_arrow_property(
                    part_id, part_type, self.specs['velocity_arrows_type_scale'],
                    self.specs['velocity_arrows_type_colors'],
                    self.specs['velocity_arrows_type_radii'], 'velocity')
                reset_material = True

            if self.specs['force_arrows']:
                self._draw_arrow_property(
                    part_id, part_type, self.specs['force_arrows_type_scale'],
                    self.specs['force_arrows_type_colors'],
                    self.specs['force_arrows_type_radii'], 'force')
                reset_material = True

            if self.specs['director_arrows']:
                self._draw_arrow_property(
                    part_id, part_type, self.specs['director_arrows_type_scale'],
                    self.specs['director_arrows_type_colors'],
                    self.specs['director_arrows_type_radii'], 'director')
                reset_material = True

    def _draw_arrow_property(self, part_id, part_type,
                             type_scale, type_colors, type_radii, prop):
        sc = self._modulo_indexing(type_scale, part_type)
        if sc > 0:
            v = self.particles[prop][part_id]
            col = self._modulo_indexing(type_colors, part_type)
            radius = self._modulo_indexing(type_radii, part_type)
            draw_arrow(
                self.particles['pos'][part_id], np.array(v, dtype=float) * sc,
                radius, col, self.materials['chrome'], self.specs['quality_arrows'])

    def _draw_bonds(self):
        box_l_2 = self.system.box_l / 2.0
        for b in self.bonds:
            col = self._modulo_indexing(self.specs['bond_type_colors'], b[2])
            mat = self.materials[self._modulo_indexing(
                self.specs['bond_type_materials'], b[2])]
            radius = self._modulo_indexing(
                self.specs['bond_type_radius'], b[2])
            x_a = self.particles['pos'][b[0]]
            x_b = self.particles['pos'][b[1]]
            dx = x_b - x_a

            if abs(dx[0]) < box_l_2[0] and abs(
                    dx[1]) < box_l_2[1] and abs(dx[2]) < box_l_2[2]:
                # BOND COMPLETELY INSIDE BOX
                draw_cylinder(x_a, x_b, radius, col, mat,
                              self.specs['quality_bonds'])
                if self.has_images:
                    for imx in self.image_vectors[0]:
                        for imy in self.image_vectors[1]:
                            for imz in self.image_vectors[2]:
                                if imx != 0 or imy != 0 or imz != 0:
                                    offset = [imx, imy, imz] * \
                                        self.system.box_l
                                    draw_cylinder(x_a + offset, x_b + offset,
                                                  radius, col, mat, self.specs['quality_bonds'])
            else:
                # BOND CROSSES THE BOX BOUNDARIES
                d = self._cut_bond(x_a, dx)
                if d is np.inf:
                    continue

                s_a = x_a + d * dx
                s_b = x_b - (1 - d) * dx
                draw_cylinder(x_a, s_a, radius, col, mat,
                              self.specs['quality_bonds'])
                draw_cylinder(x_b, s_b, radius, col, mat,
                              self.specs['quality_bonds'])

                if self.has_images:
                    for imx in self.image_vectors[0]:
                        for imy in self.image_vectors[1]:
                            for imz in self.image_vectors[2]:
                                if imx != 0 or imy != 0 or imz != 0:
                                    offset = [imx, imy, imz] * \
                                        self.system.box_l
                                    draw_cylinder(x_a + offset, s_a + offset, radius, col, mat,
                                                  self.specs['quality_bonds'])
                                    draw_cylinder(x_b + offset, s_b + offset, radius, col, mat,
                                                  self.specs['quality_bonds'])

    def _redraw_sphere(self, pos):
        OpenGL.GL.glPushMatrix()
        OpenGL.GL.glTranslatef(pos[0], pos[1], pos[2])
        OpenGL.GL.glCallList(self.dl_sphere)
        OpenGL.GL.glPopMatrix()

    # HELPER TO DRAW PERIODIC BONDS
    def _cut_bond(self, x_a, dx):
        """
        Algorithm for the line-plane intersection problem. Given the unfolded
        positions of two particles, determine: 1) the box image that minimizes
        the folded bond length, 2) which side of the simulation box is crossed
        by the bond, 3) how much of the bond is inside the unit cell starting
        from the first particle. That scalar is returned by the function, and
        can be used to calculate the coordinates of the line-plane
        intersection point. The algorithm can be found at
        https://en.wikipedia.org/w/index.php?title=Line%E2%80%93plane_intersection&oldid=940427752#Algebraic_form
        """
        # corner case: the two particles occupy the same position
        if np.dot(dx, dx) < 1e-9:
            return np.inf
        # find the box image that minimizes the unfolded bond length
        shift = np.rint(dx / self.system.box_l)
        # get the unfolded bond vector
        dx -= shift * self.system.box_l
        # find which side of the simulation box is crossed by the bond
        best_d = np.inf
        for i in range(3):
            if dx[i] == 0:
                continue  # corner case: bond is parallel to a face
            elif dx[i] > 0:
                p0_i = self.system.box_l[i]
            else:
                p0_i = 0
            # calculate the length of the bond that is inside the box using
            # an optimized version of `np.dot(p0 - x0, n) / np.dot(dx, n)`
            # where the dot products decay to an array access, because `n`
            # is always a unit vector in rectangular boxes and its sign
            # is canceled out by the division
            d = (p0_i - x_a[i]) / dx[i]
            if d < best_d:
                best_d = d
        return best_d

    # ARROWS IN A PLANE FOR LB VELOCITIES
    def _draw_lb_vel(self):

        for lbl in self.lb_plane_vel:
            p = lbl[0]
            v = lbl[1]
            draw_arrow(
                p, v *
                self.specs['LB_vel_scale'],
                self.lb_arrow_radius,
                self.specs['LB_arrow_color'],
                self.materials[self.specs['LB_arrow_material']],
                self.specs['LB_arrow_quality'])

    # USE MODULO IF THERE ARE MORE PARTICLE TYPES THAN TYPE DEFINITIONS FOR
    # COLORS, MATERIALS ETC.
    @staticmethod
    def _modulo_indexing(l, t):
        return l[t % len(l)]

    # FADE PARTICLE CHARGE COLOR FROM WHITE (q=0) to PLUSCOLOR (q=q_max) RESP
    # MINUSCOLOR (q=q_min)
    def _color_by_charge(self, q):
        if q < 0:
            c = 1.0 * q / self.min_q
            return np.array(
                self.specs['particle_charge_colors'][0]) * c + (1 - c) * np.array([1, 1, 1])
        else:
            c = 1.0 * q / self.max_q

            return np.array(
                self.specs['particle_charge_colors'][1]) * c + (1 - c) * np.array([1, 1, 1])

    # ON INITIALIZATION, CHECK q_max/q_min
    def _update_charge_color_range(self):
        if len(self.particles['charge'][:]) > 0:
            self.min_q = min(self.particles['charge'][:])
            self.max_q = max(self.particles['charge'][:])

    def _handle_screenshot(self):
        if self.take_screenshot:
            self.take_screenshot = False
            data = OpenGL.GL.glReadPixels(0, 0, self.specs['window_size'][0],
                                          self.specs['window_size'][1],
                                          OpenGL.GL.GL_RGB, OpenGL.GL.GL_FLOAT)
            script_name = os.path.splitext(sys.argv[0])[0]

            i = 0
            while os.path.exists(f"{script_name}_{i:04d}.png"):
                i += 1
            file_name = f"{script_name}_{i:04d}.png"

            data = np.flipud(data.reshape((data.shape[1], data.shape[0], 3)))
            imsave(file_name, data)

            self.screenshot_captured = True
            self.screenshot_capture_time = time.time()
            self.screenshot_capture_txt = f"Saved screenshot {file_name}"

    def _display_all(self):

        OpenGL.GL.glClear(
            OpenGL.GL.GL_COLOR_BUFFER_BIT | OpenGL.GL.GL_DEPTH_BUFFER_BIT)

        OpenGL.GL.glLoadMatrixf(self.camera.modelview)

        self._set_camera_spotlight()

        self._draw_system()
        self._draw_texts()

        OpenGL.GLUT.glutSwapBuffers()

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

            self._draw_text(10, 10, f"{self.fps} fps", self.text_color)
            self._draw_text(10, 30, f"{1000.0 / self.fps} ms/frame",
                            self.text_color)
            self.fps_count += 1

        # DRAW PARTICLE INFO
        if self.show_system_info:
            self._draw_sysinfo_dict(self.system_info)
        elif self.info_id != -1:
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
                self._draw_text(
                    self.specs['window_size'][0] - len(
                        self.screenshot_capture_txt) * 9.0 - 15,
                    self.specs['window_size'][1] - 15, self.screenshot_capture_txt, col)

    def _draw_sysinfo_dict(self, sysinfo_dict):
        y = 0
        for key, value_list in sysinfo_dict.items():
            # CATEGORY TITLE
            self._draw_text(10, self.specs['window_size'][1] - 10 - 15 * y,
                            f"{key}:", self.text_color)
            # ITEM LIST
            for item in value_list:
                txt = str(item)
                n_lines = int(
                    len(txt) * 9.0 / self.specs['window_size'][0]) + 1
                chars_per_line = int(self.specs['window_size'][0] / 9) - 4
                if chars_per_line < 20:
                    break
                ls = 0
                for _ in range(n_lines):
                    y += 1
                    line_txt = txt[ls:ls + chars_per_line]
                    self._draw_text(
                        30, self.specs['window_size'][1] - 10 - 15 * y,
                        line_txt, self.text_color)
                    ls += chars_per_line
            y += 1.5

    def _draw_particle_dict(self, particle_dict, max_len):
        y = 0
        for k, v in particle_dict.items():
            txt = f"{k} {(max_len - len(k)) * ' '} {v}"
            self._draw_text(10, self.specs['window_size'][1] - 10 - 15 * y,
                            txt, self.text_color)
            y += 1

    # CALLED ON WINDOW POSITION/SIZE CHANGE
    def _reshape_window(self, w, h):
        OpenGL.GL.glViewport(0, 0, w, h)
        OpenGL.GL.glMatrixMode(OpenGL.GL.GL_PROJECTION)
        OpenGL.GL.glLoadIdentity()
        box_diag = np.linalg.norm(self.system.box_l)
        OpenGL.GLU.gluPerspective(40, 1.0 * w / h,
                                  self.specs['close_cut_distance'],
                                  self.specs['far_cut_distance'] * box_diag)
        OpenGL.GL.glMatrixMode(OpenGL.GL.GL_MODELVIEW)
        self.specs['window_size'][0] = w
        self.specs['window_size'][1] = h

    # INITS FOR GLUT FUNCTIONS
    def _init_OpenGL_callbacks(self):
        # OpenGL Callbacks
        def display():
            if self.hasParticleData and self.glut_main_loop_started:
                self._display_all()
            return

        # pylint: disable=unused-argument
        def keyboard_up(button, x, y):
            if isinstance(button, bytes):
                button = button.decode("utf-8")
            self.keyboard_manager.keyboard_up(button)
            return

        # pylint: disable=unused-argument
        def keyboard_down(button, x, y):
            if isinstance(button, bytes):
                button = button.decode("utf-8")
            self.keyboard_manager.keyboard_down(button)
            return

        def mouse(button, state, x, y):
            self.mouse_manager.mouse_click(button, state, x, y)
            return

        def motion(x, y):
            self.mouse_manager.mouse_move(x, y)
            return

        def redraw_on_idle():
            # DON'T REPOST FASTER THAN 60 FPS
            self.draw_elapsed += (time.time() - self.draw_timer)
            if self.draw_elapsed > 1.0 / 60.0:
                self.draw_elapsed = 0
                OpenGL.GLUT.glutPostRedisplay()
            self.draw_timer = time.time()
            return

        def reshape_callback(w, h):
            self._reshape_window(w, h)

        def close_window():
            os._exit(1)

        OpenGL.GLUT.glutDisplayFunc(display)
        OpenGL.GLUT.glutMouseFunc(mouse)
        OpenGL.GLUT.glutKeyboardFunc(keyboard_down)
        OpenGL.GLUT.glutKeyboardUpFunc(keyboard_up)
        OpenGL.GLUT.glutSpecialFunc(keyboard_down)
        OpenGL.GLUT.glutSpecialUpFunc(keyboard_up)
        OpenGL.GLUT.glutReshapeFunc(reshape_callback)
        OpenGL.GLUT.glutMotionFunc(motion)
        OpenGL.GLUT.glutWMCloseFunc(close_window)

        OpenGL.GLUT.glutIdleFunc(redraw_on_idle)

    def _init_timers(self):

        # TIMERS FOR register_callback
        def dummy_timer(index):
            self.timers[index][1]()
            OpenGL.GLUT.glutTimerFunc(
                self.timers[index][0],
                dummy_timer,
                index)

        index = 0
        for t in self.timers:
            OpenGL.GLUT.glutTimerFunc(t[0], dummy_timer, index)
            index += 1

        # HANDLE INPUT WITH 60FPS
        # pylint: disable=unused-argument
        def timed_handle_input(data):
            self.keyboard_manager.handle_input()
            OpenGL.GLUT.glutTimerFunc(17, timed_handle_input, -1)

        OpenGL.GLUT.glutTimerFunc(17, timed_handle_input, -1)

    # CLICKED ON PARTICLE: DRAG; CLICKED ON BACKGROUND: CAMERA
    def _mouse_motion(self, mouse_pos, mouse_pos_old, mouse_button_state):

        if self.specs['drag_enabled'] and self.drag_id != -1:
            part_pos = self.particles['pos'][self.drag_id]
            viewport = OpenGL.GL.glGetIntegerv(OpenGL.GL.GL_VIEWPORT)
            mouse_world = OpenGL.GLU.gluUnProject(
                mouse_pos[0], viewport[3] - mouse_pos[1], self.depth)

            self.dragExtForce = self.specs['drag_force'] * \
                (np.asarray(mouse_world) - np.array(part_pos))
            self.trigger_set_particle_drag = True
        else:
            self.camera.rotate_camera(mouse_pos, mouse_pos_old,
                                      mouse_button_state)

    # DRAW SCENE AGAIN WITHOUT LIGHT TO IDENTIFY PARTICLE ID BY PIXEL COLOR
    def _get_particle_id(self, pos):

        OpenGL.GL.glClearColor(0.0, 0.0, 0.0, 1.0)
        OpenGL.GL.glClear(
            OpenGL.GL.GL_COLOR_BUFFER_BIT | OpenGL.GL.GL_DEPTH_BUFFER_BIT)

        OpenGL.GL.glLoadMatrixf(self.camera.modelview)

        OpenGL.GL.glDisable(OpenGL.GL.GL_LIGHTING)
        OpenGL.GL.glDisable(OpenGL.GL.GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            OpenGL.GL.glDisable(OpenGL.GL.GL_LIGHT1)
        self._draw_system_particles(color_by_id=True)
        viewport = OpenGL.GL.glGetIntegerv(OpenGL.GL.GL_VIEWPORT)

        read_pixel = OpenGL.GL.glReadPixelsui(
            pos[0], viewport[3] - pos[1], 1, 1,
            OpenGL.GL.GL_RGB, OpenGL.GL.GL_FLOAT)[0][0]
        depth = OpenGL.GL.glReadPixelsf(
            pos[0], viewport[3] - pos[1], 1, 1,
            OpenGL.GL.GL_DEPTH_COMPONENT, OpenGL.GL.GL_FLOAT)[0][0]

        part_id = self._fcolor_to_id(read_pixel)

        OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHTING)
        OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHT1)
        OpenGL.GL.glClearColor(self.specs['background_color'][0],
                               self.specs['background_color'][1],
                               self.specs['background_color'][2], 1.)

        return part_id, depth

    @staticmethod
    def _id_to_fcolor(part_id):
        part_id += 1
        return [int(part_id / 256 ** 2) / 255.0,
                int((part_id % 256 ** 2) / 256) / 255.0,
                (part_id % 256) / 255.0, 1.0]

    @staticmethod
    def _fcolor_to_id(fcolor):
        if (fcolor == [0, 0, 0]).all():
            return -1
        else:
            return int(fcolor[0] * 255) * 256 ** 2 + \
                int(fcolor[1] * 255) * 256 + \
                int(fcolor[2] * 255) - 1

    # pylint: disable=unused-argument
    def _set_particle_drag(self, pos, pos_old):
        part_id, depth = self._get_particle_id(pos)
        self.drag_id = part_id

        if part_id != -1:
            self.dragPosInitial = self.particles['pos'][self.drag_id]
            self.extForceOld = self.particles['ext_force'][self.drag_id][:]
            self.depth = depth

    # pylint: disable=unused-argument
    def _reset_particle_drag(self, pos, pos_old):
        if self.drag_id != -1:
            self.trigger_reset_particle_drag = True

    def _get_particle_info(self, pos):
        part_id, _ = self._get_particle_id(pos)
        if self.show_system_info:
            self.show_system_info = False
        elif part_id == -1 and self.info_id == -1:
            self.show_system_info = True
            self.update_system_info()
        self.info_id = part_id

    def _next_particle_info(self):
        self.info_id = (self.info_id + 1) % len(self.particles['pos'])

    def _previous_particle_info(self):
        self.info_id = (self.info_id - 1) % len(self.particles['pos'])

    # ESPRESSO RELATED INITS
    def _init_espresso_visualization(self):
        self.max_q = 0
        self.min_q = 0

        self.drag_id = -1
        self.info_id = -1
        self.show_system_info = False
        self.dragPosInitial = []
        self.extForceOld = []
        self.dragExtForceOld = []
        self.trigger_reset_particle_drag = False
        self.trigger_set_particle_drag = False

        self.depth = 0

        # LOOK FOR LB ACTOR
        if self.specs['LB_draw_velocity_plane'] or \
                self.specs['LB_draw_nodes'] or \
                self.specs['LB_draw_node_boundaries']:
            lb_types = [espressomd.lb.LBFluid]
            if espressomd.has_features('CUDA'):
                lb_types.append(espressomd.lb.LBFluidGPU)
            for a in self.system.actors:
                if isinstance(a, tuple(lb_types)):
                    # if 'agrid' in pa:
                    self.lb_params = a.get_params()
                    self.lb = a
                    self.lb_is_cpu = isinstance(a, espressomd.lb.LBFluid)
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

        self.box_eqn = np.full((6, 4), np.nan)
        # set face normals of box
        self.box_eqn[:3, :3] = np.identity(3)
        self.box_eqn[3:, :3] = -np.identity(3)
        # set face parameters of box
        self.box_eqn[:3, 3] = 0.001 * self.system.box_l
        self.box_eqn[3:, 3] = 1.001 * self.system.box_l

        self.camera.center = np.array(self.system.box_l) * 0.5
        self.camera.update_modelview()

    # DEFAULT CONTROLS
    def _init_controls(self):

        # MOUSE LOOK/ROTATE/DRAG
        self.mouse_manager.register_button(MouseButtonEvent(
            None, MouseFireEvent.FreeMotion, self._mouse_motion))

        self.mouse_manager.register_button(MouseButtonEvent(
            3, MouseFireEvent.ButtonPressed, self.camera.move_backward))

        self.mouse_manager.register_button(MouseButtonEvent(
            4, MouseFireEvent.ButtonPressed, self.camera.move_forward))

        # START/STOP DRAG
        if self.specs['drag_enabled']:
            self.mouse_manager.register_button(MouseButtonEvent(
                OpenGL.GLUT.GLUT_LEFT_BUTTON, MouseFireEvent.ButtonPressed, self._set_particle_drag, True))
            self.mouse_manager.register_button(MouseButtonEvent(
                OpenGL.GLUT.GLUT_LEFT_BUTTON, MouseFireEvent.ButtonReleased, self._reset_particle_drag, True))

        # PARTICLE INFORMATION
        self.mouse_manager.register_button(MouseButtonEvent(
            OpenGL.GLUT.GLUT_LEFT_BUTTON, MouseFireEvent.DoubleClick, self._get_particle_info, True))

        # CYCLE THROUGH PARTICLES
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            OpenGL.GLUT.GLUT_KEY_LEFT, KeyboardFireEvent.Pressed, self._previous_particle_info))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            OpenGL.GLUT.GLUT_KEY_RIGHT, KeyboardFireEvent.Pressed, self._next_particle_info))

        # <SPACE> PAUSE INTEGRATION THREAD
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            ' ', KeyboardFireEvent.Pressed, self._pause, True))

        # <RETURN> TAKE SCREENSHOT
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            '\x0d', KeyboardFireEvent.Pressed, self._trigger_screenshot, True))

        # <ESCAPE> QUIT
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            '\x1b', KeyboardFireEvent.Pressed, self._quit, True))

        # CAMERA CONTROL VIA KEYBOARD
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'w', KeyboardFireEvent.Hold, self.camera.move_up, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            's', KeyboardFireEvent.Hold, self.camera.move_down, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'a', KeyboardFireEvent.Hold, self.camera.move_left, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'd', KeyboardFireEvent.Hold, self.camera.move_right, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'e', KeyboardFireEvent.Hold, self.camera.rotate_system_XR, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'q', KeyboardFireEvent.Hold, self.camera.rotate_system_XL, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'c', KeyboardFireEvent.Hold, self.camera.rotate_system_YR, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'z', KeyboardFireEvent.Hold, self.camera.rotate_system_YL, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'r', KeyboardFireEvent.Hold, self.camera.rotate_system_ZR, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'f', KeyboardFireEvent.Hold, self.camera.rotate_system_ZL, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            't', KeyboardFireEvent.Hold, self.camera.move_forward, True))
        self.keyboard_manager.register_button(KeyboardButtonEvent(
            'g', KeyboardFireEvent.Hold, self.camera.move_backward, True))

    # CALLED ON ESCAPE PRESSED. TRIGGERS sys.exit() after ES is done
    def _quit(self):
        self.quit_safely = True

    def _pause(self):
        self.paused = not self.paused

    def _trigger_screenshot(self):
        self.take_screenshot = True

    # ASYNCHRONOUS PARALLEL CALLS OF glLight CAUSES SEG FAULTS, SO ONLY CHANGE
    # LIGHT AT CENTRAL display METHOD AND TRIGGER CHANGES
    def _set_camera_spotlight(self):
        if self.specs['spotlight_enabled']:
            p = self.camera.cam_pos
            fp = [p[0], p[1], p[2], 1]
            OpenGL.GL.glLightfv(OpenGL.GL.GL_LIGHT1, OpenGL.GL.GL_POSITION, fp)
            OpenGL.GL.glLightfv(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_SPOT_DIRECTION,
                self.camera.state_target)

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

        self._set_camera_spotlight()

        return Camera(cam_pos=np.array(cp), cam_target=ct, cam_right=cr,
                      move_speed=0.5 * box_diag / 17.0, center=box_center)

    def _init_opengl(self):
        OpenGL.GLUT.glutInit(self.specs['name'])
        OpenGL.GLUT.glutInitDisplayMode(
            OpenGL.GLUT.GLUT_DOUBLE | OpenGL.GLUT.GLUT_RGB | OpenGL.GLUT.GLUT_DEPTH)

        OpenGL.GLUT.glutInitWindowSize(self.specs['window_size'][0],
                                       self.specs['window_size'][1])

        OpenGL.GLUT.glutCreateWindow(b"ESPResSo visualization")

        OpenGL.GL.glClearColor(self.specs['background_color'][0],
                               self.specs['background_color'][1],
                               self.specs['background_color'][2], 1.)

        OpenGL.GL.glEnable(OpenGL.GL.GL_BLEND)
        OpenGL.GL.glBlendFunc(
            OpenGL.GL.GL_SRC_ALPHA,
            OpenGL.GL.GL_ONE_MINUS_SRC_ALPHA)

        OpenGL.GL.glEnable(OpenGL.GL.GL_POINT_SMOOTH)
        OpenGL.GL.glEnable(OpenGL.GL.GL_LINE_SMOOTH)
        OpenGL.GL.glHint(OpenGL.GL.GL_LINE_SMOOTH_HINT, OpenGL.GL.GL_NICEST)

        # BAD FOR TRANSPARENT PARTICLES
        # OpenGL.GL.glEnable(OpenGL.GL.GL_CULL_FACE)
        # OpenGL.GL.glCullFace(OpenGL.GL.GL_BACK)

        OpenGL.GL.glLineWidth(2.0)
        OpenGL.GLUT.glutIgnoreKeyRepeat(1)

        # setup lighting
        if self.specs['light_size'] == 'auto':
            box_diag = np.linalg.norm(self.system.box_l)
            self.specs['light_size'] = box_diag * 2.0

        OpenGL.GL.glEnable(OpenGL.GL.GL_DEPTH_TEST)
        OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHTING)

        # LIGHT0
        if self.specs['light_pos'] != 'auto':
            OpenGL.GL.glLightfv(OpenGL.GL.GL_LIGHT0, OpenGL.GL.GL_POSITION,
                                np.array(self.specs['light_pos']).tolist())
        else:
            OpenGL.GL.glLightfv(OpenGL.GL.GL_LIGHT0, OpenGL.GL.GL_POSITION,
                                (np.array(self.system.box_l) * 1.1).tolist())

        OpenGL.GL.glLightfv(
            OpenGL.GL.GL_LIGHT0,
            OpenGL.GL.GL_AMBIENT,
            self.specs['light_colors'][0])
        OpenGL.GL.glLightfv(
            OpenGL.GL.GL_LIGHT0,
            OpenGL.GL.GL_DIFFUSE,
            self.specs['light_colors'][1])
        OpenGL.GL.glLightfv(
            OpenGL.GL.GL_LIGHT0,
            OpenGL.GL.GL_SPECULAR,
            self.specs['light_colors'][2])

        OpenGL.GL.glLightf(
            OpenGL.GL.GL_LIGHT0, OpenGL.GL.GL_CONSTANT_ATTENUATION,
            1.0 / self.specs['light_brightness'])
        OpenGL.GL.glLightf(
            OpenGL.GL.GL_LIGHT0, OpenGL.GL.GL_LINEAR_ATTENUATION,
            1.0 / self.specs['light_size'])
        OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHT0)

        # LIGHT1: SPOTLIGHT ON CAMERA IN LOOK DIRECTION
        if self.specs['spotlight_enabled']:
            OpenGL.GL.glLightfv(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_POSITION,
                [0, 0, 0, 1])

            OpenGL.GL.glLightfv(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_AMBIENT,
                self.specs['spotlight_colors'][0])
            OpenGL.GL.glLightfv(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_DIFFUSE,
                self.specs['spotlight_colors'][1])
            OpenGL.GL.glLightfv(OpenGL.GL.GL_LIGHT1, OpenGL.GL.GL_SPECULAR,
                                self.specs['spotlight_colors'][2])

            OpenGL.GL.glLightf(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_SPOT_CUTOFF,
                self.specs['spotlight_angle'])
            OpenGL.GL.glLightfv(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_SPOT_DIRECTION,
                [1.0, 1.0, 1.0])
            OpenGL.GL.glLightf(OpenGL.GL.GL_LIGHT1, OpenGL.GL.GL_SPOT_EXPONENT,
                               self.specs['spotlight_focus'])

            OpenGL.GL.glLightf(
                OpenGL.GL.GL_LIGHT1, OpenGL.GL.GL_CONSTANT_ATTENUATION,
                1.0 / self.specs['spotlight_brightness'])
            OpenGL.GL.glLightf(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_LINEAR_ATTENUATION,
                0.0)
            OpenGL.GL.glLightf(
                OpenGL.GL.GL_LIGHT1,
                OpenGL.GL.GL_QUADRATIC_ATTENUATION,
                0.0)
            OpenGL.GL.glEnable(OpenGL.GL.GL_LIGHT1)

        self.dl_sphere = OpenGL.GL.glGenLists(1)


# END OF MAIN CLASS

class Shape():
    """
    Shape base class in the visualizer context.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        self.shape = shape
        self.particle_type = particle_type
        self.color = color
        self.material = material
        self.quality = quality
        self.box_l = box_l
        self.rasterize_resolution = rasterize_resolution
        self.pointsize = rasterize_pointsize

        self.rasterized_surface_points = None

    def draw(self):
        """
        Draw shape via rasterization. Used as a default draw method.
        Can and should be overwritten in child classes to implement a better
        draw method.

        """
        # get and store points of rasterized surface if not already present
        if self.rasterized_surface_points is None:
            self.rasterized_surface_points = self._rasterize_shape()

        set_solid_material(self.color, self.material)
        OpenGL.GL.glPointSize(self.pointsize)
        OpenGL.GL.glBegin(OpenGL.GL.GL_POINTS)
        for point in self.rasterized_surface_points:
            OpenGL.GL.glVertex3f(point[0], point[1], point[2])
        OpenGL.GL.glEnd()

    def _rasterize_shape(self):
        # rasterize brute force
        spacing = max(self.box_l) / self.rasterize_resolution
        resolution = np.array(self.box_l) / spacing

        points = []
        for i in range(int(resolution[0])):
            for j in range(int(resolution[1])):
                for k in range(int(resolution[2])):
                    # some shapes may not have a well-defined distance function in the whole domain
                    # and may throw upon asking for a distance
                    try:
                        p = np.array([i, j, k]) * spacing
                        dist, vec = self.shape.call_method(
                            "calc_distance", position=p.tolist())
                        if not np.isnan(vec).any() and not np.isnan(
                                dist) and abs(dist) < spacing:
                            points.append((p - vec).tolist())
                    # domain error translates to ValueError (cython)
                    except ValueError:
                        continue
        return points


class Cylinder(Shape):
    """
    Drawable Shape Cylinder.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = np.array(self.shape.get_parameter('center'))
        self.axis = np.array(self.shape.get_parameter('axis'))
        self.length = self.shape.get_parameter('length')
        self.radius = self.shape.get_parameter('radius')
        self.cap_center_1 = self.center - self.axis / \
            np.linalg.norm(self.axis) * 0.5 * self.length
        self.cap_center_2 = self.center + self.axis / \
            np.linalg.norm(self.axis) * 0.5 * self.length

    def draw(self):
        draw_cylinder(self.cap_center_1, self.cap_center_2,
                      self.radius, self.color, self.material,
                      self.quality, draw_caps=True)


class Ellipsoid(Shape):
    """
    Drawable Shape Ellipsoid.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = np.array(self.shape.get_parameter('center'))
        self.semiaxis_a = np.array(self.shape.get_parameter('a'))
        self.semiaxis_b = np.array(self.shape.get_parameter('b'))
        self.semiaxis_c = np.array(self.shape.get_parameter('b'))

    def draw(self):
        set_solid_material(self.color, self.material)
        OpenGL.GL.glPushMatrix()
        OpenGL.GL.glTranslatef(self.center[0], self.center[1], self.center[2])
        OpenGL.GL.glScalef(self.semiaxis_a, self.semiaxis_b, self.semiaxis_c)
        OpenGL.GLUT.glutSolidSphere(1, self.quality, self.quality)
        OpenGL.GL.glPopMatrix()


class HollowConicalFrustum(Shape):
    """
    Drawable Shape HollowConicalFrustum.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = np.array(self.shape.get_parameter('center'))
        self.radius_1 = np.array(self.shape.get_parameter('r1'))
        self.radius_2 = np.array(self.shape.get_parameter('r2'))
        self.length = np.array(self.shape.get_parameter('length'))
        self.thickness = np.array(self.shape.get_parameter('thickness'))
        self.axis = np.array(self.shape.get_parameter('axis'))

    def draw(self):
        """
        Draw using OpenGL Extrusion library, if available.
        Use rasterization of base class, otherwise.

        """
        if bool(OpenGL.GLE.gleSpiral):
            self._draw_using_gle()
        else:
            super().draw()

    # if available, use the GL Extrusion library
    def _draw_using_gle(self):
        set_solid_material(self.color, self.material)
        OpenGL.GL.glPushMatrix()
        # basic position and orientation
        OpenGL.GL.glTranslate(self.center[0], self.center[1], self.center[2])
        ax, rx, ry = rotation_helper(self.axis)
        OpenGL.GL.glRotatef(ax, rx, ry, 0.0)

        n = max(20, self.quality)
        rotation_angle = np.arctan(
            (self.radius_1 - self.radius_2) / self.length)
        l = self.length / np.cos(rotation_angle)

        contour = []
        for theta in np.linspace(-0.5 * np.pi, 0.5 * np.pi, n):
            contour.append([np.sin(theta) * 0.5 * self.thickness,
                            np.cos(theta) * 0.5 * self.thickness + 0.5 * l])
        for theta in np.linspace(0.5 * np.pi, -0.5 * np.pi, n):
            contour.append([np.sin(theta) * 0.5 * self.thickness,
                            -np.cos(theta) * 0.5 * self.thickness - 0.5 * l])
        contour = np.matmul(np.array(contour),
                            np.array([[np.cos(rotation_angle), -np.sin(rotation_angle)],
                                      [np.sin(rotation_angle), np.cos(rotation_angle)]]))

        normals = np.diff(np.array(contour), axis=0)
        normals /= np.linalg.norm(normals, ord=2, axis=1, keepdims=True)
        normals = np.roll(normals, 1, axis=1)

        OpenGL.GLE.gleSetJoinStyle(OpenGL.GLE.TUBE_JN_ANGLE)
        OpenGL.GLE.gleSetNumSides(max(90, 3 * self.quality))
        OpenGL.GLE.gleSpiral(contour, normals, [0, 0, 1], 0.5 * (self.radius_1 + self.radius_2), 0., 0., 0.,
                             [[1, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0]], 0., 360)

        OpenGL.GL.glPopMatrix()


class SimplePore(Shape):
    """
    Drawable Shape SimplePore.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = np.array(self.shape.get_parameter('center'))
        self.axis = np.array(self.shape.get_parameter('axis'))
        self.length = np.array(self.shape.get_parameter('length'))
        self.radius = np.array(self.shape.get_parameter('radius'))
        self.smoothing_radius = np.array(
            self.shape.get_parameter('smoothing_radius'))
        self.max_box_l = max(box_l)

    def draw(self):
        """
        Draw using OpenGL Extrusion library, if available.
        Use OpenGL primitives + clip planes, otherwise.

        """
        set_solid_material(self.color, self.material)
        OpenGL.GL.glPushMatrix()
        # basic position and orientation
        OpenGL.GL.glTranslate(self.center[0], self.center[1], self.center[2])
        ax, rx, ry = rotation_helper(self.axis)
        OpenGL.GL.glRotatef(ax, rx, ry, 0.0)

        if bool(OpenGL.GLE.gleSpiral):
            self._draw_using_gle()
        else:
            self._draw_using_primitives()

        OpenGL.GL.glPopMatrix()

    def _draw_using_gle(self):
        # if available, use the GL Extrusion library
        if bool(OpenGL.GLE.gleSpiral):
            n = max(10, self.quality // 3)
            contour = [[0.5 * self.max_box_l, -0.5 * self.length]]
            for theta in np.linspace(0, 0.5 * np.pi, n):
                contour.append([(1. - np.sin(theta)) * self.smoothing_radius,
                                -0.5 * self.length + (1. - np.cos(theta)) * self.smoothing_radius])
            for theta in np.linspace(0.5 * np.pi, np.pi, n):
                contour.append([(1. - np.sin(theta)) * self.smoothing_radius,
                                0.5 * self.length - (1. + np.cos(theta)) * self.smoothing_radius])
            contour.append([0.5 * self.max_box_l, 0.5 * self.length])

            normals = np.diff(np.array(contour), axis=0)
            normals /= np.linalg.norm(normals, ord=2, axis=1, keepdims=True)
            normals = np.roll(normals, 1, axis=1)
            normals[:, 0] *= -1

            OpenGL.GLE.gleSetJoinStyle(OpenGL.GLE.TUBE_JN_ANGLE)
            OpenGL.GLE.gleSetNumSides(max(90, 3 * self.quality))
            OpenGL.GLE.gleSpiral(contour, normals, [0, 0, 1], self.radius, 0., 0., 0.,
                                 [[1, 0, 0], [0, 1, 0]], [[0, 0, 0], [0, 0, 0]], 0., 360)

    def _draw_using_primitives(self):
        clip_plane = get_extra_clip_plane()
        # cylinder
        OpenGL.GL.glTranslate(0, 0, -0.5 * self.length + self.smoothing_radius)
        quadric = OpenGL.GLU.gluNewQuadric()
        OpenGL.GLU.gluCylinder(quadric, self.radius, self.radius, self.length - 2 *
                               self.smoothing_radius, self.quality, self.quality)
        # torus segment
        OpenGL.GL.glEnable(clip_plane)
        OpenGL.GL.glClipPlane(clip_plane, (0, 0, -1, 0))
        OpenGL.GLUT.glutSolidTorus(
            self.smoothing_radius,
            self.radius + self.smoothing_radius,
            self.quality,
            self.quality)
        OpenGL.GL.glDisable(clip_plane)
        # wall
        OpenGL.GL.glTranslate(0, 0, -self.smoothing_radius)
        OpenGL.GLU.gluPartialDisk(quadric, self.radius + self.smoothing_radius,
                                  2.0 * self.max_box_l, self.quality, 1, 0, 360)
        # torus segment
        OpenGL.GL.glTranslate(0, 0, self.length - self.smoothing_radius)
        OpenGL.GL.glEnable(clip_plane)
        OpenGL.GL.glClipPlane(clip_plane, (0, 0, 1, 0))
        OpenGL.GLUT.glutSolidTorus(
            self.smoothing_radius,
            self.radius + self.smoothing_radius,
            self.quality,
            self.quality)
        OpenGL.GL.glDisable(clip_plane)
        # wall
        OpenGL.GL.glTranslate(0, 0, self.smoothing_radius)
        OpenGL.GLU.gluPartialDisk(quadric, self.radius + self.smoothing_radius,
                                  2.0 * self.max_box_l, self.quality, 1, 0, 360)

        OpenGL.GLU.gluDeleteQuadric(quadric)


class Slitpore(Shape):
    """
    Drawable Shape Slitpore.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.channel_width = np.array(
            self.shape.get_parameter('channel_width'))
        self.lower_smoothing_radius = np.array(
            self.shape.get_parameter('lower_smoothing_radius'))
        self.upper_smoothing_radius = np.array(
            self.shape.get_parameter('upper_smoothing_radius'))
        self.pore_length = np.array(self.shape.get_parameter('pore_length'))
        self.pore_mouth = np.array(self.shape.get_parameter('pore_mouth'))
        self.pore_width = np.array(self.shape.get_parameter('pore_width'))
        self.max_box_l = max(self.box_l)

    def draw(self):
        set_solid_material(self.color, self.material)
        # If pore is large, an additional wall is necessary
        if self.pore_width > 2. * self.lower_smoothing_radius:
            wall_0 = [
                [0.5 * (self.max_box_l - self.pore_width) + self.lower_smoothing_radius,
                 0., self.pore_mouth - self.pore_length],
                [0.5 * (self.max_box_l + self.pore_width) - self.lower_smoothing_radius,
                 0., self.pore_mouth - self.pore_length],
                [0.5 * (self.max_box_l + self.pore_width) - self.lower_smoothing_radius,
                 self.max_self.box_l, self.pore_mouth - self.pore_length],
                [0.5 * (self.max_box_l - self.pore_width) + self.lower_smoothing_radius, self.max_box_l, self.pore_mouth - self.pore_length]]
            draw_plane(wall_0, self.color, self.material)

        # Add the remaining walls
        wall_1 = [
            [0., 0., self.channel_width + self.pore_mouth],
            [self.max_box_l, 0., self.channel_width + self.pore_mouth],
            [self.max_box_l, self.max_box_l, self.channel_width + self.pore_mouth],
            [0., self.max_box_l, self.channel_width + self.pore_mouth]]

        wall_2 = [
            [0., 0., self.pore_mouth],
            [0.5 * (self.max_box_l - self.pore_width) -
             self.upper_smoothing_radius, 0., self.pore_mouth],
            [0.5 * (self.max_box_l - self.pore_width) -
             self.upper_smoothing_radius, self.max_box_l, self.pore_mouth],
            [0., self.max_box_l, self.pore_mouth]]

        wall_3 = [
            [0.5 * (self.max_box_l + self.pore_width) +
             self.upper_smoothing_radius, 0., self.pore_mouth],
            [self.max_box_l, 0., self.pore_mouth],
            [self.max_box_l, self.max_box_l, self.pore_mouth],
            [0.5 * (self.max_box_l + self.pore_width) + self.upper_smoothing_radius, self.max_box_l, self.pore_mouth]]

        wall_4 = [
            [0.5 * (self.max_box_l - self.pore_width), 0.,
             self.pore_mouth - self.upper_smoothing_radius],
            [0.5 * (self.max_box_l - self.pore_width), self.max_box_l,
             self.pore_mouth - self.upper_smoothing_radius],
            [0.5 * (self.max_box_l - self.pore_width), self.max_box_l, self.pore_mouth -
             self.pore_length + self.lower_smoothing_radius],
            [0.5 * (self.max_box_l - self.pore_width), 0., self.pore_mouth - self.pore_length + self.lower_smoothing_radius]]

        wall_5 = [
            [0.5 * (self.max_box_l + self.pore_width), 0.,
             self.pore_mouth - self.upper_smoothing_radius],
            [0.5 * (self.max_box_l + self.pore_width), self.max_box_l,
             self.pore_mouth - self.upper_smoothing_radius],
            [0.5 * (self.max_box_l + self.pore_width), self.max_box_l, self.pore_mouth -
             self.pore_length + self.lower_smoothing_radius],
            [0.5 * (self.max_box_l + self.pore_width), 0., self.pore_mouth - self.pore_length + self.lower_smoothing_radius]]

        draw_plane(wall_1, self.color, self.material)
        draw_plane(wall_2, self.color, self.material)
        draw_plane(wall_3, self.color, self.material)
        draw_plane(wall_4, self.color, self.material)
        draw_plane(wall_5, self.color, self.material)

        # Add smooth edges via clipped cylinders
        ax, rx, ry = rotation_helper([0., 1., 0.])

        OpenGL.GL.glPushMatrix()
        quadric = OpenGL.GLU.gluNewQuadric()
        OpenGL.GL.glTranslate(0.5 * self.max_box_l - self.upper_smoothing_radius -
                              0.5 * self.pore_width, 0, self.pore_mouth - self.upper_smoothing_radius)
        OpenGL.GL.glRotatef(ax, rx, ry, 0.)

        # Upper edges
        clip_plane = get_extra_clip_plane()
        OpenGL.GL.glEnable(clip_plane)
        OpenGL.GL.glClipPlane(clip_plane,
                              (1, -1, 0, -self.upper_smoothing_radius))
        OpenGL.GLU.gluCylinder(quadric, self.upper_smoothing_radius,
                               self.upper_smoothing_radius, self.max_box_l, self.quality, self.quality)

        OpenGL.GL.glTranslate(self.pore_width + 2. * self.upper_smoothing_radius,
                              0, 0)
        OpenGL.GL.glClipPlane(clip_plane,
                              (-1, -1, 0, -self.upper_smoothing_radius))
        OpenGL.GLU.gluCylinder(quadric, self.upper_smoothing_radius,
                               self.upper_smoothing_radius, self.max_box_l, self.quality, self.quality)

        # Lower edges
        OpenGL.GL.glTranslate(- self.upper_smoothing_radius - self.lower_smoothing_radius,
                              self.pore_length - self.upper_smoothing_radius - self.lower_smoothing_radius, 0)
        OpenGL.GL.glClipPlane(clip_plane,
                              (1, 1, 0, -self.lower_smoothing_radius))
        OpenGL.GLU.gluCylinder(quadric, self.lower_smoothing_radius,
                               self.lower_smoothing_radius, self.max_box_l, self.quality, self.quality)

        OpenGL.GL.glTranslate(-self.pore_width + 2. *
                              self.lower_smoothing_radius, 0, 0)
        OpenGL.GL.glClipPlane(clip_plane,
                              (-1, 1, 0, -self.lower_smoothing_radius))
        OpenGL.GLU.gluCylinder(quadric, self.lower_smoothing_radius,
                               self.lower_smoothing_radius, self.max_box_l, self.quality, self.quality)

        OpenGL.GL.glDisable(clip_plane)
        OpenGL.GLU.gluDeleteQuadric(quadric)
        OpenGL.GL.glPopMatrix()


class Sphere(Shape):
    """
    Drawable Shape Sphere.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = self.shape.get_parameter('center')
        self.radius = self.shape.get_parameter('radius')

    def draw(self):
        OpenGL.GL.glPushMatrix()
        OpenGL.GL.glTranslatef(self.center[0], self.center[1], self.center[2])
        set_solid_material(self.color, self.material)
        OpenGL.GLUT.glutSolidSphere(self.radius, self.quality, self.quality)
        OpenGL.GL.glPopMatrix()


class Spherocylinder(Shape):
    """
    Drawable Shape Spherocylinder.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.center = np.array(self.shape.get_parameter('center'))
        self.axis = np.array(self.shape.get_parameter('axis'))
        self.length = self.shape.get_parameter('length')
        self.radius = self.shape.get_parameter('radius')
        self.cap_center_1 = self.center - self.axis / \
            np.linalg.norm(self.axis) * 0.5 * self.length
        self.cap_center_2 = self.center + self.axis / \
            np.linalg.norm(self.axis) * 0.5 * self.length

    def draw(self):
        set_solid_material(self.color, self.material)
        OpenGL.GL.glPushMatrix()
        quadric = OpenGL.GLU.gluNewQuadric()

        d = self.cap_center_2 - self.cap_center_1
        if d[2] == 0.0:
            d[2] = 0.0001

        v = np.linalg.norm(d)
        if v == 0:
            ax = 57.2957795
        else:
            ax = 57.2957795 * math.acos(d[2] / v)

        if d[2] < 0.0:
            ax = -ax
        rx = -d[1] * d[2]
        ry = d[0] * d[2]
        length = np.linalg.norm(d)
        OpenGL.GL.glTranslatef(
            self.cap_center_1[0],
            self.cap_center_1[1],
            self.cap_center_1[2])
        OpenGL.GL.glRotatef(ax, rx, ry, 0.0)

        # First hemispherical cap
        clip_plane = get_extra_clip_plane()
        OpenGL.GL.glEnable(clip_plane)
        OpenGL.GL.glClipPlane(clip_plane, (0, 0, -1, 0))
        OpenGL.GLU.gluSphere(quadric, self.radius, self.quality, self.quality)
        OpenGL.GL.glDisable(clip_plane)
        # Cylinder
        OpenGL.GLU.gluCylinder(quadric, self.radius, self.radius,
                               length, self.quality, self.quality)
        # Second hemispherical cap
        OpenGL.GL.glTranslatef(0, 0, v)
        OpenGL.GL.glEnable(clip_plane)
        OpenGL.GL.glClipPlane(clip_plane, (0, 0, 1, 0))
        OpenGL.GLU.gluSphere(quadric, self.radius, self.quality, self.quality)
        OpenGL.GL.glDisable(clip_plane)

        OpenGL.GLU.gluDeleteQuadric(quadric)
        OpenGL.GL.glPopMatrix()


class Wall(Shape):
    """
    Drawable Shape Wall.

    """

    def __init__(self, shape, particle_type, color, material,
                 quality, box_l, rasterize_resolution, rasterize_pointsize):
        super().__init__(shape, particle_type, color, material,
                         quality, box_l, rasterize_resolution, rasterize_pointsize)
        self.distance = self.shape.get_parameter('dist')
        self.normal = self.shape.get_parameter('normal')
        self.box_diag = np.linalg.norm(self.box_l)
        self.edges = self._edges_from_pn(self.distance * np.array(self.normal),
                                         self.normal, 2 * self.box_diag)

    @staticmethod
    def _get_tangents(n):
        n = np.array(n)
        v1 = np.random.randn(3)
        v1 -= v1.dot(n) * n / np.linalg.norm(n)**2
        v2 = np.cross(n, v1)
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        return v1, v2

    def _edges_from_pn(self, p, n, diag):
        v1, v2 = self._get_tangents(n)
        edges = [p + diag * v1, p + diag * v2, p - diag * v1, p - diag * v2]
        return edges

    def draw(self):
        draw_plane(self.edges, self.color, self.material)


# OPENGL DRAW WRAPPERS

def set_solid_material(color, material=(0.6, 1.0, 0.1, 0.4, 1.0)):
    OpenGL.GL.glMaterialfv(OpenGL.GL.GL_FRONT_AND_BACK, OpenGL.GL.GL_AMBIENT, [
        color[0] * material[0], color[1] * material[0], color[2] * material[0], material[4]])
    OpenGL.GL.glMaterialfv(OpenGL.GL.GL_FRONT_AND_BACK, OpenGL.GL.GL_DIFFUSE, [
        color[0] * material[1], color[1] * material[1], color[2] * material[1], material[4]])
    OpenGL.GL.glMaterialfv(OpenGL.GL.GL_FRONT_AND_BACK, OpenGL.GL.GL_SPECULAR, [
        color[0] * material[2], color[1] * material[2], color[2] * material[2], material[4]])
    OpenGL.GL.glMaterialf(
        OpenGL.GL.GL_FRONT_AND_BACK,
        OpenGL.GL.GL_SHININESS,
        int(material[3] * 128))


def draw_box(p0, s, color, material, width):
    OpenGL.GL.glLineWidth(width)
    set_solid_material(color, material)
    OpenGL.GL.glPushMatrix()
    OpenGL.GL.glTranslatef(p0[0], p0[1], p0[2])
    OpenGL.GL.glBegin(OpenGL.GL.GL_LINE_LOOP)
    OpenGL.GL.glVertex3f(0.0, 0.0, 0.0)
    OpenGL.GL.glVertex3f(s[0], 0.0, 0.0)
    OpenGL.GL.glVertex3f(s[0], s[1], 0.0)
    OpenGL.GL.glVertex3f(0, s[1], 0.0)
    OpenGL.GL.glEnd()
    OpenGL.GL.glBegin(OpenGL.GL.GL_LINE_LOOP)
    OpenGL.GL.glVertex3f(0.0, 0.0, s[2])
    OpenGL.GL.glVertex3f(s[0], 0.0, s[2])
    OpenGL.GL.glVertex3f(s[0], s[1], s[2])
    OpenGL.GL.glVertex3f(0, s[1], s[2])
    OpenGL.GL.glEnd()
    OpenGL.GL.glBegin(OpenGL.GL.GL_LINES)
    OpenGL.GL.glVertex3f(0.0, 0.0, 0.0)
    OpenGL.GL.glVertex3f(0.0, 0.0, s[2])
    OpenGL.GL.glVertex3f(s[0], 0.0, 0.0)
    OpenGL.GL.glVertex3f(s[0], 0.0, s[2])
    OpenGL.GL.glVertex3f(s[0], s[1], 0.0)
    OpenGL.GL.glVertex3f(s[0], s[1], s[2])
    OpenGL.GL.glVertex3f(0.0, s[1], 0.0)
    OpenGL.GL.glVertex3f(0.0, s[1], s[2])
    OpenGL.GL.glEnd()

    OpenGL.GL.glPopMatrix()


def draw_plane(corners, color, material):

    set_solid_material(color, material)

    OpenGL.GL.glBegin(OpenGL.GL.GL_QUADS)
    for c in corners:
        OpenGL.GL.glVertex3f(c[0], c[1], c[2])
    OpenGL.GL.glEnd()


def draw_cylinder(posA, posB, radius, color, material, quality,
                  draw_caps=False):
    set_solid_material(color, material)
    OpenGL.GL.glPushMatrix()
    quadric = OpenGL.GLU.gluNewQuadric()

    d = posB - posA

    # angle,t,length = calcAngle(d)
    length = np.linalg.norm(d)
    OpenGL.GL.glTranslatef(posA[0], posA[1], posA[2])

    ax, rx, ry = rotation_helper(d)
    OpenGL.GL.glRotatef(ax, rx, ry, 0.0)
    OpenGL.GLU.gluCylinder(quadric, radius, radius, length, quality, quality)

    if draw_caps:
        OpenGL.GLU.gluDisk(quadric, 0, radius, quality, quality)
        OpenGL.GL.glTranslatef(0, 0, length)
        OpenGL.GLU.gluDisk(quadric, 0, radius, quality, quality)

    OpenGL.GLU.gluDeleteQuadric(quadric)
    OpenGL.GL.glPopMatrix()


def rotation_helper(d):
    # the rotation axis is the cross product between z and d
    vz = np.cross([0.0, 0.0, 1.0], d)
    # get the angle using a dot product
    angle = 180.0 / np.pi * math.acos(d[2] / np.linalg.norm(d))

    return angle, vz[0], vz[1]


def get_extra_clip_plane():

    # ON SOME HARDWARE (e.g. MAC) only 6 CLIP PLANES ARE ALLOWED,
    # SO CAPPING OF BOX BOUNDARIES AND ADDITIONAL SHAPE CLIP PLANES
    # ARE NOT POSSIBLE. THIS WILL CAUSE THE SHAPES THAT NEED ADDITIONAL
    # CLIP PLANES TO NOT BE CLIPPED ON ONE FACE OF THE BOX

    if sys.platform == "darwin":
        return OpenGL.GL.GL_CLIP_PLANE0
    else:
        return OpenGL.GL.GL_CLIP_PLANE0 + 6


def draw_arrow(pos, d, radius, color, material, quality):
    pos2 = np.array(pos) + np.array(d)

    draw_cylinder(pos, pos2, radius, color, material, quality)

    ax, rx, ry = rotation_helper(d)

    OpenGL.GL.glPushMatrix()
    OpenGL.GL.glTranslatef(pos2[0], pos2[1], pos2[2])
    OpenGL.GL.glRotatef(ax, rx, ry, 0.0)
    OpenGL.GLUT.glutSolidCone(radius * 3, radius * 3, quality, quality)
    OpenGL.GL.glPopMatrix()


# MOUSE EVENT MANAGER
class MouseFireEvent:
    """Event type of mouse button used for mouse callbacks.

    """

    ButtonPressed = 0
    FreeMotion = 1
    ButtonMotion = 2
    ButtonReleased = 3
    DoubleClick = 4


class MouseButtonEvent:
    """Mouse event used for mouse callbacks. Stores button and callback.

    """

    def __init__(self, button, fireEvent, callback, positional=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.positional = positional


class MouseManager:
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
        self.mouseState = {OpenGL.GLUT.GLUT_LEFT_BUTTON: OpenGL.GLUT.GLUT_UP,
                           OpenGL.GLUT.GLUT_MIDDLE_BUTTON: OpenGL.GLUT.GLUT_UP,
                           OpenGL.GLUT.GLUT_RIGHT_BUTTON: OpenGL.GLUT.GLUT_UP,
                           '3': OpenGL.GLUT.GLUT_UP, '4': OpenGL.GLUT.GLUT_UP}
        self.pressedTime = {OpenGL.GLUT.GLUT_LEFT_BUTTON: 0,
                            OpenGL.GLUT.GLUT_MIDDLE_BUTTON: 0,
                            OpenGL.GLUT.GLUT_RIGHT_BUTTON: 0,
                            3: 0, 4: 0}
        self.pressedTimeOld = {OpenGL.GLUT.GLUT_LEFT_BUTTON: 0, OpenGL.GLUT.GLUT_MIDDLE_BUTTON: 0,
                               OpenGL.GLUT.GLUT_RIGHT_BUTTON: 0, 3: 0, 4: 0}

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
            if me.button == button and state == OpenGL.GLUT.GLUT_DOWN:
                if me.positional:
                    me.callback(self.mousePos, self.mousePosOld)
                else:
                    me.callback()
        for me in self.mouseEventsReleased:
            if me.button == button and state == OpenGL.GLUT.GLUT_UP:
                if me.positional:
                    me.callback(self.mousePos, self.mousePosOld)
                else:
                    me.callback()

        if state == OpenGL.GLUT.GLUT_DOWN:
            self.pressedTimeOld[button] = self.pressedTime[button]
            self.pressedTime[button] = time.time()

        for me in self.mouseEventsDoubleClick:
            if me.button == button and state == OpenGL.GLUT.GLUT_DOWN and self.pressedTime[
                    button] - self.pressedTimeOld[button] < 0.25:
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

class KeyboardFireEvent:
    """Event type of button used for keyboard callbacks.

    """

    Pressed = 0
    Hold = 1
    Released = 2


class KeyboardButtonEvent:
    """Keyboard event used for keyboard callbacks. Stores button, event type and callback.

    """

    def __init__(self, button, fireEvent, callback, internal=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.internal = internal


class KeyboardManager:
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
        if button not in self.keyStateOld.keys():
            self.keyStateOld[button] = 0


# CAMERA

class Camera:

    def __init__(self, cam_pos=np.array([0, 0, 1]), cam_target=np.array([0, 0, 0]),
                 cam_right=np.array([1.0, 0.0, 0.0]), move_speed=0.5,
                 global_rot_speed=3.0, center=np.array([0, 0, 0])):
        self.cam_pos = cam_pos
        self.move_speed = move_speed
        self.global_rot_speed = global_rot_speed
        self.center = center

        self.modelview = np.identity(4, np.float32)

        t = cam_pos - cam_target
        r = np.linalg.norm(t)

        self.state_target = -t / r
        self.state_right = cam_right / np.linalg.norm(cam_right)
        self.state_up = np.cross(self.state_right, self.state_target)
        self.state_pos = np.array([0, 0, r])

        self.update_modelview()

    def move_forward(self):
        self.state_pos[2] += self.move_speed
        self.update_modelview()

    def move_backward(self):
        self.state_pos[2] -= self.move_speed
        self.update_modelview()

    def move_up(self):
        self.state_pos[1] += self.move_speed
        self.update_modelview()

    def move_down(self):
        self.state_pos[1] -= self.move_speed
        self.update_modelview()

    def move_left(self):
        self.state_pos[0] -= self.move_speed
        self.update_modelview()

    def move_right(self):
        self.state_pos[0] += self.move_speed
        self.update_modelview()

    def rotate_system_XL(self):
        self.rotate_system_y(0.01 * self.global_rot_speed)

    def rotate_system_XR(self):
        self.rotate_system_y(-0.01 * self.global_rot_speed)

    def rotate_system_YL(self):
        self.rotate_system_z(0.01 * self.global_rot_speed)

    def rotate_system_YR(self):
        self.rotate_system_z(-0.01 * self.global_rot_speed)

    def rotate_system_ZL(self):
        self.rotate_system_x(0.01 * self.global_rot_speed)

    def rotate_system_ZR(self):
        self.rotate_system_x(-0.01 * self.global_rot_speed)

    def rotate_camera(self, mouse_pos, mouse_pos_old, mouse_button_state):
        dm = mouse_pos - mouse_pos_old

        if mouse_button_state[OpenGL.GLUT.GLUT_LEFT_BUTTON] == OpenGL.GLUT.GLUT_DOWN:
            if dm[0] != 0:
                self.rotate_system_y(-0.001 * dm[0] * self.global_rot_speed)
            if dm[1] != 0:
                self.rotate_system_x(-0.001 * dm[1] * self.global_rot_speed)
        elif mouse_button_state[OpenGL.GLUT.GLUT_RIGHT_BUTTON] == OpenGL.GLUT.GLUT_DOWN:
            self.state_pos[0] -= 0.05 * dm[0] * self.move_speed
            self.state_pos[1] += 0.05 * dm[1] * self.move_speed
            self.update_modelview()
        elif mouse_button_state[OpenGL.GLUT.GLUT_MIDDLE_BUTTON] == OpenGL.GLUT.GLUT_DOWN:
            self.state_pos[2] += 0.05 * dm[1] * self.move_speed
            self.rotate_system_z(dm[0] * 0.001 * self.global_rot_speed)

    def get_camera_rotation_matrix(self, target_vec, up_vec):
        normalized_target_vec = target_vec / np.linalg.norm(target_vec)
        normalized_up_vec = up_vec / np.linalg.norm(up_vec)
        perpendicular_normal = np.cross(normalized_up_vec, target_vec)

        rotation_matrix = np.identity(4, np.float32)
        rotation_matrix[:3, :3] = np.array(
            [perpendicular_normal,
             np.cross(normalized_target_vec, perpendicular_normal),
             normalized_target_vec]).T

        return rotation_matrix

    def rotate_vector(self, vector, phi, axis):
        """
        Rotate vector around (unit vector) axis by angle phi.
        Uses Rodrigues' rotation formula.

        """
        rotated = vector * np.cos(phi) + np.cross(axis, vector) * np.sin(phi) \
            + axis * np.dot(axis, vector) * (1 - np.cos(phi))
        return rotated

    def rotate_system_x(self, angle):
        self.state_target = self.rotate_vector(self.state_target,
                                               angle, self.state_right)
        self.state_up = np.cross(self.state_right, self.state_target)
        self.update_modelview()

    def rotate_system_y(self, angle):
        self.state_target = self.rotate_vector(self.state_target,
                                               angle, self.state_up)
        self.state_right = np.cross(self.state_target, self.state_up)
        self.update_modelview()

    def rotate_system_z(self, angle):
        self.state_right = self.rotate_vector(self.state_right,
                                              angle, self.state_target)
        self.state_up = self.rotate_vector(self.state_up,
                                           angle, self.state_target)
        self.update_modelview()

    def update_modelview(self):

        self.state_up /= np.linalg.norm(self.state_up)
        self.state_right /= np.linalg.norm(self.state_right)
        self.state_target /= np.linalg.norm(self.state_target)

        # Center Box
        trans = np.identity(4, np.float32)
        trans[3, :3] -= self.center

        # Camera rotation
        rotate_cam = self.get_camera_rotation_matrix(
            -self.state_target, self.state_up)

        # System translation
        trans_cam = np.identity(4, np.float32)
        trans_cam[3, :3] -= self.state_pos

        self.modelview = trans.dot(rotate_cam.dot(trans_cam))

        c_xyz = -1 * np.mat(self.modelview[:3, :3]) * \
            np.mat(self.modelview[3, :3]).T
        self.cam_pos = np.array([c_xyz[0, 0], c_xyz[1, 0], c_xyz[2, 0]])
