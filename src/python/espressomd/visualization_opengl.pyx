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

import scipy.spatial
include "myconfig.pxi"


class openGLLive(object):

    def __init__(self, system, **kwargs):

        # DEFAULT PROPERTIES
        self.specs = {
            'window_size':				  [800, 800],
            'name':						  'Espresso Visualization',
            'background_color':			  [0, 0, 0],
            'update_fps':				  30.0,
            'periodic_images':	   		  [0, 0, 0],
            'draw_box':		 	  		  True,
            'quality_particles':          16,
            'quality_bonds':              16,
            'quality_arrows':             16,
            'quality_constraints':        32,
            'close_cut_distance':         0.1,
            'far_cut_distance':           5,
            'camera_position':	   		  'auto',
            'camera_rotation':	   		  [3.55, -0.4],

            'particle_coloring':   		  'auto',
            'particle_sizes':			  'auto',
            'particle_type_colors':		  [[1, 1, 0, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            'particle_type_materials':	  [[0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1]],
            'particle_charge_colors':	  [[1, 0, 0, 1], [0, 1, 0, 1]],

            'draw_constraints':			  True,
            'rasterize_pointsize':	      10,
            'rasterize_resolution':	      75.0,
            'constraint_type_colors':     [[0.5, 0.5, 0.5, 0.9], [0, 0.5, 0.5, 0.9], [0.5, 0, 0.5, 0.9], [0.5, 0.5, 0, 0.9], [0, 0, 0.5, 0.9], [0.5, 0, 0, 0.9]],
            'constraint_type_materials':  [[0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1]],

            'draw_bonds':				  True,
            'bond_type_radius':		      [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            'bond_type_colors':			  [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            'bond_type_materials':	      [[0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1]],

            'ext_force_arrows': 		  True,
            'ext_force_arrows_scale': 	  [1, 1, 1, 1, 1, 1, 1],

            'LB':						  False,

            'light_pos':				  'auto',
            'light_colors':				  [[0.1, 0.1, 0.2, 1.0], [0.9, 0.9, 0.9, 1.0], [1.0, 1.0, 1.0, 1.0]],
            'light_brightness':   		  1.0,
            'light_size':		   		  'auto',

            'spotlight_enabled':          True,
            'spotlight_colors':		      [[0.1, 0.1, 0.2, 1.0], [0.5, 0.5, 0.5, 1.0], [1.0, 1.0, 1.0, 1.0]],
            'spotlight_angle':            45,
            'spotlight_focus':            1,
            'spotlight_brightness':   	  0.4,

            'drag_enabled':		   		  False,
            'drag_force':		   		  3.0
        }

        # OVERWRITE WITH USER PROPERTIES
        for key in kwargs.iterkeys():
            if key not in self.specs.iterkeys():
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
        self.keyboardManager = KeyboardManager()
        self.mouseManager = MouseManager()
        self.timers = []

    # CALLBACKS FOR THE MAIN THREAD
    def registerCallback(self, cb, interval=1000):
        self.timers.append((int(interval), cb))

    # THE BLOCKING START METHOD
    def start(self):
        self.initOpenGL()
        self.initEspressoVisualization()
        self.initCamera()
        self.initControls()
        self.initCallbacks()

        # POST DISPLAY WITH 60FPS
        def timed_update_redraw(data):
            glutPostRedisplay()
            self.keyboardManager.handleInput()
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

    # CALLED FROM ESPRESSO INTEGRATION LOOP
    # CHANGES OF ESPRESSO SYSTEM CAN ONLY HAPPEN HERE
    def update(self):
        if self.started:

            # UPDATE ON STARTUP
            if not self.hasParticleData:
                self.updateParticles()
                self.updateChargeColorRange()
                self.updateBonds()
                IF CONSTRAINTS:
                    self.updateConstraints()
                self.hasParticleData = True

            # IF CALLED TOO OFTEN, ONLY UPDATE WITH GIVEN FREQ
            self.elapsedTime += (time.time() - self.measureTimeBeforeIntegrate)
            if self.elapsedTime > 1.0 / self.specs['update_fps']:
                self.elapsedTime = 0
                self.updateParticles()
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

    # GET THE PARTICLE DATA
    def updateParticles(self):
        IF EXTERNAL_FORCES and ELECTROSTATICS:
            self.particles = {'coords':  	self.system.part[:].pos_folded,
                              'types':   	self.system.part[:].type,
                              'ext_forces': self.system.part[:].ext_force,
                              'charges':    self.system.part[:].q}
        ELIF EXTERNAL_FORCES and not ELECTROSTATICS:
            self.particles = {'coords':  	self.system.part[:].pos_folded,
                              'types':   	self.system.part[:].type,
                              'ext_forces': self.system.part[:].ext_force,
                              'charges':    [0] * len(self.system.part)}
        ELIF not EXTERNAL_FORCES and ELECTROSTATICS:
            self.particles = {'coords':  	self.system.part[:].pos_folded,
                              'types':   	self.system.part[:].type,
                              'ext_forces': [0, 0, 0] * len(self.system.part),
                              'charges':    self.system.part[:].q}
        ELIF not EXTERNAL_FORCES and not ELECTROSTATICS:
            self.particles = {'coords':  	self.system.part[:].pos_folded,
                              'types':   	self.system.part[:].type,
                              'ext_forces': [0, 0, 0] * len(self.system.part),
                              'charges':    [0] * len(self.system.part)}

    def edgesFromPN(self, p, n, diag):
        v1, v2 = self.getTangents(n)

        edges = []
        edges.append(p + diag * v1)
        edges.append(p + diag * v2)
        edges.append(p - diag * v1)
        edges.append(p - diag * v2)
        return edges

    # GET THE CONSTRAINT DATA
    def updateConstraints(self):

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
            edges = self.edgesFromPN(d * np.array(n), n, 2 * box_diag)
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
                [self.rasterizeBruteForce(s[0]), s[1]])

    def getTangents(self, n):
        if n[0] > 0.5 or n[1] > 0.5:
            v1 = np.array([n[1], -n[0], 0])
        else:
            v1 = np.array([-n[2], 0, n[0]])
        v2 = np.cross(n, v1)
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)

        return v1, v2

    def rasterizeBruteForce(self, shape):
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
    def updateBonds(self):
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
            self.drawLBVel()
        if self.specs['draw_box']:
            self.drawSystemBox()

        self.drawSystemParticles()

        if self.specs['draw_bonds']:
            self.drawBonds()

        IF CONSTRAINTS:
            if self.specs['draw_constraints']:
                self.drawConstraints()

    def drawSystemBox(self):
        drawBox([0, 0, 0], self.system.box_l, self.invBackgroundCol)

    def drawConstraints(self):

        for i in range(6):
            glEnable(GL_CLIP_PLANE0 + i)
            glClipPlane(GL_CLIP_PLANE0 + i, self.box_eqn[i])

        for s in self.shapes['Shapes::Sphere']:
            drawSphere(s[0], s[1], self.modulo_indexing(self.specs['constraint_type_colors'], s[2]), self.modulo_indexing(
                self.specs['constraint_type_materials'], s[2]), self.specs['quality_constraints'])

        for s in self.shapes['Shapes::SpheroCylinder']:
            drawSpheroCylinder(
                s[0], s[1], s[2], self.modulo_indexing(
                    self.specs['constraint_type_colors'], s[3]),
                               self.modulo_indexing(self.specs['constraint_type_materials'], s[3]), self.specs['quality_constraints'])

        for s in self.shapes['Shapes::Wall']:
            drawPlane(
                s[0], self.modulo_indexing(
                    self.specs['constraint_type_colors'], s[1]),
                      self.modulo_indexing(self.specs['constraint_type_materials'], s[1]))

        for s in self.shapes['Shapes::Cylinder']:
            drawCylinder(s[0], s[1], s[2], self.modulo_indexing(self.specs['constraint_type_colors'], s[3]), self.modulo_indexing(
                self.specs['constraint_type_materials'], s[3]), self.specs['quality_constraints'], True)

        for s in self.shapes['Shapes::Misc']:
            drawPoints(s[0], self.specs['rasterize_pointsize'],  self.modulo_indexing(
                self.specs['constraint_type_colors'], s[1]), self.modulo_indexing(self.specs['constraint_type_materials'], s[1]))

        for i in range(6):
            glDisable(GL_CLIP_PLANE0 + i)

    def determine_radius(self, ptype):
        def radiusByLJ(ptype):
            try:
                radius = self.system.non_bonded_inter[ptype, ptype].lennard_jones.get_params()[
                    'sigma'] * 0.5
            except:
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
            except:
                radius = self.radiusByLJ(ptype)
        return radius

    def drawSystemParticles(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        for pid in pIds:
            pos = coords[pid]
            q = self.particles['charges'][pid]
            ptype = self.particles['types'][pid]
            ext_f = self.particles['ext_forces'][pid]

            radius = self.determine_radius(ptype)

            material = self.modulo_indexing(
                self.specs['particle_type_materials'], ptype)

            if self.specs['particle_coloring'] == 'id':
                color = self.IdToColorf(pid)
                glColor(color)
            elif self.specs['particle_coloring'] == 'auto':
                # Color auto: Charge then Type
                if q != 0:
                    color = self.colorByCharge(q)
                else:
                    color = self.modulo_indexing(
                        self.specs['particle_type_colors'], ptype)
            elif self.specs['particle_coloring'] == 'charge':
                color = self.colorByCharge(q)
            elif self.specs['particle_coloring'] == 'type':
                color = self.modulo_indexing(
                    self.specs['particle_type_colors'], ptype)

            drawSphere(pos, radius, color, material,
                       self.specs['quality_particles'])
            for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                    for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                        if imx != 0 or imy != 0 or imz != 0:
                            redrawSphere(
                                pos + (imx * self.imPos[0] + imy * self.imPos[1] + imz * self.imPos[2]), radius, self.specs['quality_particles'])

            if self.specs['ext_force_arrows'] or pid == self.dragId:
                if ext_f[0] != 0 or ext_f[1] != 0 or ext_f[2] != 0:
                    if pid == self.dragId:
                        sc = 1
                    else:
                        sc = self.modulo_indexing(
                            self.specs['ext_force_arrows_scale'], ptype)
                    if sc > 0:
                        drawArrow(pos, np.array(ext_f) * sc, 0.25 *
                                  sc, [1, 1, 1], self.specs['quality_arrows'])

    def drawBonds(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        b2 = self.system.box_l[0] / 2.0
        box_l2_sqr = pow(b2, 2.0)
        for b in self.bonds:
            col = self.modulo_indexing(self.specs['bond_type_colors'], b[2])
            mat = self.modulo_indexing(self.specs['bond_type_materials'], b[2])
            radius = self.modulo_indexing(self.specs['bond_type_radius'], b[2])
            d = coords[b[0]] - coords[b[1]]
            bondLen_sqr = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]

            if bondLen_sqr < box_l2_sqr:
                drawCylinder(coords[b[0]], coords[b[1]], radius,
                             col, mat, self.specs['quality_bonds'])
                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                im = np.array([imx, imy, imz])
                                drawCylinder(coords[b[0]] + im * self.imPos[dim], coords[b[1]] +
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
                    if self.isInsideBox(s):
                        if lineBoxNDot < 0:
                            s0 = s
                        else:
                            s1 = s
                        hits += 1
                        if hits >= 2:
                            break
                drawCylinder(coords[b[0]], s0, radius, col,
                             mat, self.specs['quality_bonds'])
                drawCylinder(coords[b[1]], s1, radius, col,
                             mat, self.specs['quality_bonds'])

                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                im = np.array([imx, imy, imz])
                                drawCylinder(coords[b[0]] + im * self.imPos[dim], s0 + im *
                                             self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])
                                drawCylinder(coords[b[1]] + im * self.imPos[dim], s1 + im *
                                             self.imPos[dim], radius, col, mat, self.specs['quality_bonds'])

    # HELPER TO DRAW PERIODIC BONDS
    def isInsideBox(self, p):
        eps = 1e-5
        for i in range(3):
            if p[i] < -eps or p[i] > eps + self.system.box_l[i]:
                return False
        return True

    # VOXELS FOR LB VELOCITIES
    def drawLBVel(self):
        grid = 10
        velRelax = 0.2
        cubeSize = grid * 0.25
        r = np.array([grid] * 3)

        min_vel_new = np.array([1e100] * 3)
        max_vel_new = np.array([-1e100] * 3)
        for ix in range(r[0]):
            for iy in range(r[1]):
                for iz in range(r[2]):
                    c = self.system.box_l * \
                        (np.array([ix, iy, iz]) +
                         np.array([0.5, 0.5, 0.5])) / r
                    v = self.system.actors[0].lbnode_get_node_velocity(c)
                    col = (np.array(v) - self.lb_min_vel) / self.lb_vel_range
                    alpha = 0.1  # np.linalg.norm(col)
                    drawCube(c, cubeSize, col, alpha)

    # USE MODULO IF THERE ARE MORE PARTICLE TYPES THAN TYPE DEFINITIONS FOR
    # COLORS, MATERIALS ETC..
    def modulo_indexing(self, l, t):
        return l[t % len(l)]

    # FADE PARTICE CHARGE COLOR FROM WHITE (q=0) to PLUSCOLOR (q=q_max) RESP
    # MINUSCOLOR (q=q_min)
    def colorByCharge(self, q):
        if q < 0:
            c = 1.0 * q / self.minq
            return np.array(self.specs['particle_charge_colors'][0]) * c + (1 - c) * np.array([1, 1, 1, 1])
        else:
            c = 1.0 * q / self.maxq

            return np.array(self.specs['particle_charge_colors'][1]) * c + (1 - c) * np.array([1, 1, 1, 1])

    # ON INITIALIZATION, CHECK q_max/q_min
    def updateChargeColorRange(self):
        if len(self.particles['charges'][:]) > 0:
            self.minq = min(self.particles['charges'][:])
            self.maxq = max(self.particles['charges'][:])

    # INITS FOR GLUT FUNCTIONS
    def initCallbacks(self):
        # OpenGl Callbacks
        def display():
            if self.hasParticleData:
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
                glLoadIdentity()

                self.camera.glLookAt()

                if self.updateLightPos:
                    self.setLightPos()
                    self.updateLightPos = False

                self.camera.rotateSystem()

                self.draw()

                glutSwapBuffers()
            return

        def keyboardUp(button, x, y):
            self.keyboardManager.keyboardUp(button)
            return

        def keyboardDown(button, x, y):
            self.keyboardManager.keyboardDown(button)
            return

        def mouse(button, state, x, y):
            self.mouseManager.mouseClick(button, state, x, y)
            return

        def motion(x, y):
            self.mouseManager.mouseMove(x, y)
            return

        def idleUpdate():
            return

        # CALLED ION WINDOW POSITION/SIZE CHANGE
        def reshapeWindow(w, h):
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

        # TIMERS FOR registerCallback
        def dummyTimer(index):
            self.timers[index][1]()
            glutTimerFunc(self.timers[index][0], dummyTimer, index)

        glutDisplayFunc(display)
        glutMouseFunc(mouse)
        glutKeyboardFunc(keyboardDown)
        glutKeyboardUpFunc(keyboardUp)
        glutReshapeFunc(reshapeWindow)
        glutMotionFunc(motion)

        glutIdleFunc(idleUpdate)

        index = 0
        for t in self.timers:
            glutTimerFunc(t[0], dummyTimer, index)
            index += 1

    # CLICKED ON PARTICLE: DRAG; CLICKED ON BACKGROUND: CAMERA
    def mouseMotion(self, mousePos, mousePosOld, mouseButtonState):

        if self.dragId != -1:
            ppos = self.particles['coords'][self.dragId]
            viewport = glGetIntegerv(GL_VIEWPORT)
            mouseWorld = gluUnProject(
                mousePos[0], viewport[3] - mousePos[1], self.depth)

            self.dragExtForce = self.specs['drag_force'] * \
                (np.asarray(mouseWorld) - np.array(ppos))
            self.triggerSetParticleDrag = True
        else:
            self.camera.rotateCamera(mousePos, mousePosOld, mouseButtonState)

    # DRAW SCENE AGAIN WITHOUT LIGHT TO IDENTIFY PARTICLE ID BY PIXEL COLOR
    def setParticleDrag(self, pos, pos_old):

        glClearColor(0.0, 0.0, 0.0, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        self.camera.glLookAt()
        self.camera.rotateSystem()

        oldColMode = self.specs['particle_coloring']
        self.specs['particle_coloring'] = 'id'
        glDisable(GL_LIGHTING)
        glDisable(GL_LIGHT0)
        if self.specs['spotlight_enabled']:
            glDisable(GL_LIGHT1)
        self.drawSystemParticles()
        viewport = glGetIntegerv(GL_VIEWPORT)

        readPixel = glReadPixelsui(
            pos[0], viewport[3] - pos[1], 1, 1, GL_RGB, GL_FLOAT)[0][0]
        depth = glReadPixelsf(
            pos[0], viewport[3] - pos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]
        pid = self.fcolorToId(readPixel)

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

    def resetParticleDrag(self, pos, pos_old):
        if self.dragId != -1:
            self.triggerResetParticleDrag = True

    def IdToColorf(self, pid):
        pid += 1
        return [int(pid / (256 * 256)) / 255.0, int((pid % (256 * 256)) / 256) / 255.0, (pid % 256) / 255.0, 1.0]

    def fcolorToId(self, fcol):
        if (fcol == [0, 0, 0]).all():
            return -1
        else:
            return 256 * 256 * int(fcol[0] * 255) + 256 * int(fcol[1] * 255) + int(fcol[2] * 255) - 1

    # ALL THE INITS
    def initEspressoVisualization(self):
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

        self.lb_min_vel = np.array([-1e-6] * 3)
        self.lb_max_vel = np.array([1e-6] * 3)
        self.lb_vel_range = self.lb_max_vel - self.lb_min_vel
        self.lb_min_dens = np.array([0] * 3)
        self.lb_max_dens = np.array([0] * 3)

        self.elapsedTime = 0
        self.measureTimeBeforeIntegrate = 0

        self.boxSizeDependence()

    # BOX PLANES (NORMAL, ORIGIN) FOR PERIODIC BONDS
    def boxSizeDependence(self):
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
    def initControls(self):
        # MOUSE LOOK/ROTATE/DRAG
        self.mouseManager.registerButton(MouseButtonEvent(
            None, MouseFireEvent.FreeMotion, self.mouseMotion))

        self.mouseManager.registerButton(MouseButtonEvent(
            3, MouseFireEvent.ButtonPressed, self.camera.moveForward))

        self.mouseManager.registerButton(MouseButtonEvent(
            4, MouseFireEvent.ButtonPressed, self.camera.moveBackward))

        # START/STOP DRAG
        if self.specs['drag_enabled']:
            self.mouseManager.registerButton(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonPressed, self.setParticleDrag, True))
            self.mouseManager.registerButton(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonReleased, self.resetParticleDrag, True))

        # KEYBOARD BUTTONS
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'w', KeyboardFireEvent.Hold, self.camera.moveForward, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            's', KeyboardFireEvent.Hold, self.camera.moveBackward, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'a', KeyboardFireEvent.Hold, self.camera.moveLeft, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'd', KeyboardFireEvent.Hold, self.camera.moveRight, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'e', KeyboardFireEvent.Hold, self.camera.rotateSystemXR, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'q', KeyboardFireEvent.Hold, self.camera.rotateSystemXL, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'z', KeyboardFireEvent.Hold, self.camera.rotateSystemYR, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'c', KeyboardFireEvent.Hold, self.camera.rotateSystemYL, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'r', KeyboardFireEvent.Hold, self.camera.rotateSystemZR, True))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'f', KeyboardFireEvent.Hold, self.camera.rotateSystemZL, True))

    # ASYNCHRONOUS PARALLEL CALLS OF glLight CAUSES SEG FAULTS, SO ONLY CHANGE
    # LIGHT AT CENTRAL display METHOD AND TRIGGER CHANGES
    def setLightPos(self):
        if self.specs['light_pos'] == 'auto':
            glLightfv(GL_LIGHT0, GL_POSITION, [
                      self.smooth_light_pos[0], self.smooth_light_pos[1], self.smooth_light_pos[2], 0.6])

        self.setCameraSpotlight()

    def setCameraSpotlight(self):
        p = self.camera.camPos + 10 * self.camera.lookDir
        fp = [p[0], p[1], p[2], 1]
        glLightfv(GL_LIGHT1, GL_POSITION, fp)
        glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, self.camera.lookDir)

    def triggerLightPosUpdate(self):
        self.updateLightPos = True

    def initCamera(self):
        box_diag = pow(pow(self.system.box_l[0], 2) + pow(
            self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)
        bl = self.system.box_l[0]
        bl2 = bl / 2.0
        if self.specs['camera_position'] == 'auto':
            cp = [bl * 1.3, bl * 1.3, bl * 2.5]
        else:
            cp = self.specs['camera_position']
        box_center = np.array([bl2, bl2, bl2])
        self.camera = Camera(camPos=np.array(cp), camRot=np.array(
            self.specs['camera_rotation']), moveSpeed=0.5 * box_diag / 17.0,  center=box_center, updateLights=self.triggerLightPosUpdate)
        self.smooth_light_pos = np.copy(box_center)
        self.smooth_light_posV = np.array([0.0, 0.0, 0.0])
        self.particle_COM = np.copy(box_center)
        self.setCameraSpotlight()
        self.updateLightPos = True

    def initOpenGL(self):
        glutInit(self.specs['name'])
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
        glutInitWindowSize(self.specs['window_size'][
                           0], self.specs['window_size'][1])
        glutCreateWindow(self.specs['name'])
        glClearColor(self.specs['background_color'][0], self.specs[
                     'background_color'][1], self.specs['background_color'][2], 1.)

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        glEnable(GL_BLEND)

        #glEnable(GL_CULL_FACE)

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
            glLightfv(GL_LIGHT0, GL_POSITION, self.specs['light_pos'])
        else:
            glLightfv(GL_LIGHT0, GL_POSITION, self.system.box_l * 0.5)

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


def setSolidMaterial(r, g, b, a=1.0, ambient=0.6, diffuse=1.0, specular=0.1):
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  [
                 ambient * r, ambient * g, ambient * g, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  [
                 diffuse * r, diffuse * g, diffuse * b, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  [
                 specular * r, specular * g, specular * g, a])
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50)


def setOutlineMaterial(r, g, b, a=1.0):
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  [r, g, b, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  [0, 0, 0, a])
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, [0, 0, 0, a])
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100)


def drawBox(p0, s, color):
    setSolidMaterial(color[0], color[1], color[2], 1, 2, 1)
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


def drawSphere(pos, radius, color, material, quality):
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    setSolidMaterial(color[0], color[1], color[2], color[
                     3], material[0], material[1], material[2])
    glutSolidSphere(radius, quality, quality)
    glPopMatrix()


def redrawSphere(pos, radius, quality):
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glutSolidSphere(radius, quality, quality)
    glPopMatrix()


def drawPlane(edges, color, material):

    setSolidMaterial(color[0], color[1], color[2], color[
                     3], material[0], material[1], material[2])

    glBegin(GL_QUADS)
    for e in edges:
        glVertex3f(e[0], e[1], e[2])
    glEnd()


def drawCube(pos, size, color, alpha):
    setSolidMaterial(color[0], color[1], color[2], alpha)
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glutSolidCube(size)
    glPopMatrix()


def calcAngle(d):
    z = np.array([0, 0, 1])
    l = np.linalg.norm(d)
    dot_z = np.dot(z, d)
    t = np.cross(z, d)
    return 180.0 / pi * acos(dot_z / l), t, l


def drawTriangles(triangles, color, material):
    np.random.seed(1)

    glBegin(GL_TRIANGLES)
    for t in triangles:
        color = np.random.random(3).tolist()
        color.append(1)
        setSolidMaterial(color[0], color[1], color[2],
                         color[3], material[0], material[1], material[2])
        for p in t:
            glVertex3f(p[0], p[1], p[2])
    glEnd()


def drawPoints(points, pointsize, color, material):
    setSolidMaterial(color[0], color[1], color[2], color[3],
                     material[0], material[1], material[2])
    glEnable(GL_POINT_SMOOTH)

    glPointSize(pointsize)
    glBegin(GL_POINTS)
    for p in points:
        glVertex3f(p[0], p[1], p[2])
    glEnd()


def drawCylinder(posA, posB, radius, color, material, quality, draw_caps=False):
    setSolidMaterial(color[0], color[1], color[2], color[3],
                     material[0], material[1], material[2])
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

    # angle,t,length = calcAngle(d)
    length = np.linalg.norm(d)
    glTranslatef(posA[0], posA[1], posA[2])

    glRotatef(ax, rx, ry, 0.0)
    gluCylinder(quadric, radius, radius, length, quality, quality)

    if draw_caps:
        gluDisk(quadric, 0, radius, quality, quality)
        glTranslatef(0, 0, v)
        gluDisk(quadric, 0, radius, quality, quality)

    glPopMatrix()


def drawSpheroCylinder(posA, posB, radius, color, material, quality):
    setSolidMaterial(color[0], color[1], color[
                     2], color[3], material[0], material[1], material[2])
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
    glEnable(GL_CLIP_PLANE0)
    glClipPlane(GL_CLIP_PLANE0, (0, 0, -1, 0))
    gluSphere(quadric, radius, quality, quality)
    glDisable(GL_CLIP_PLANE0)
    # Cylinder
    gluCylinder(quadric, radius, radius, length, quality, quality)
    # Second hemispherical cap
    glTranslatef(0, 0, v)
    glEnable(GL_CLIP_PLANE0)
    glClipPlane(GL_CLIP_PLANE0, (0, 0, 1, 0))
    gluSphere(quadric, radius, quality, quality)
    glDisable(GL_CLIP_PLANE0)

    glPopMatrix()


def drawArrow(pos, d, radius, color, quality):
    setSolidMaterial(color[0], color[1], color[2])
    glPushMatrix()
    glPushMatrix()

    v = np.linalg.norm(d)
    ax = 57.2957795 * acos(d[2] / v)
    if d[2] < 0.0:
        ax = -ax
    rx = -d[1] * d[2]
    ry = d[0] * d[2]

    # angle,t,length = calcAngle(d)
    glTranslatef(pos[0], pos[1], pos[2])
    glRotatef(ax, rx, ry, 0.0)
    quadric = gluNewQuadric()
    gluCylinder(quadric, radius, radius, v, quality, quality)

    glPopMatrix()
    e = pos + d
    glTranslatef(e[0], e[1], e[2])
    glRotatef(ax, rx, ry, 0.0)

    glutSolidCone(radius * 3, radius * 3, quality, quality)
    glPopMatrix()


# MOUSE EVENT MANAGER
class MouseFireEvent(object):

    ButtonPressed = 0
    FreeMotion = 1
    ButtonMotion = 2
    ButtonReleased = 3


class MouseButtonEvent(object):

    def __init__(self, button, fireEvent, callback, positional=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.positional = positional


class MouseManager(object):

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

    def registerButton(self, mouseEvent):
        if mouseEvent.fireEvent == MouseFireEvent.ButtonPressed:
            self.mouseEventsPressed.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonReleased:
            self.mouseEventsReleased.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.FreeMotion:
            self.mouseEventsFreeMotion.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonMotion:
            self.mouseEventsButtonMotion.append(mouseEvent)

    def mouseClick(self, button, state, x, y):
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

    def mouseMove(self, x, y):
        self.mousePosOld = self.mousePos
        self.mousePos = np.array([x, y])

        for me in self.mouseEventsFreeMotion:
            me.callback(self.mousePos, self.mousePosOld, self.mouseState)

# KEYBOARD EVENT MANAGER


class KeyboardFireEvent(object):

    Pressed = 0
    Hold = 1
    Released = 2


class KeyboardButtonEvent(object):

    def __init__(self, button, fireEvent, callback, internal=False):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback
        self.internal = internal


class KeyboardManager(object):

    def __init__(self):
        self.pressedKeys = set([])
        self.keyStateOld = {}
        self.keyState = {}
        self.buttonEventsPressed = []
        self.buttonEventsHold = []
        self.buttonEventsReleased = []
        self.userCallbackStack = []

    def registerButton(self, buttonEvent):
        if buttonEvent.fireEvent == KeyboardFireEvent.Pressed:
            self.buttonEventsPressed.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Hold:
            self.buttonEventsHold.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Released:
            self.buttonEventsReleased.append(buttonEvent)

    def callbackOnButton(self, be, b):
        if be.button == b:
            if be.internal:
                be.callback()
            else:
                self.userCallbackStack.append(be.callback)

    def handleInput(self):
        removeKeys = set([])
        for b in self.pressedKeys:
            if self.keyStateOld[b] == 0 and self.keyState[b] == 1:
                for be in self.buttonEventsPressed:
                    self.callbackOnButton(be, b)
                for be in self.buttonEventsHold:
                    self.callbackOnButton(be, b)

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 1:
                for be in self.buttonEventsHold:
                    self.callbackOnButton(be, b)

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 0:
                for be in self.buttonEventsReleased:
                    self.callbackOnButton(be, b)
                removeKeys.add(b)

            self.keyStateOld[b] = self.keyState[b]

        self.pressedKeys = self.pressedKeys.difference(removeKeys)

    def keyboardUp(self, button):
        self.keyState[button] = 0  # Key up

    def keyboardDown(self, button):
        self.pressedKeys.add(button)
        self.keyState[button] = 1  # Key down
        if not button in self.keyStateOld.keys():
            self.keyStateOld[button] = 0

# CAMERA


class Camera(object):

    def __init__(self, camPos=np.array([0, 0, 1]), camRot=np.array([pi, 0]), moveSpeed=0.5, rotSpeed=0.001, globalRotSpeed=3.0, center=np.array([0, 0, 0]), updateLights=None):
        self.moveSpeed = moveSpeed
        self.lookSpeed = rotSpeed
        self.globalRotSpeed = globalRotSpeed
        self.camPos = camPos
        self.camRot = camRot
        self.center = center
        self.camRotGlobal = np.array([0.0, 0.0, 0.0])
        self.updateLights = updateLights
        self.calcCameraDirections()

    def moveForward(self):
        self.camPos += self.lookDir * self.moveSpeed
        self.updateLights()

    def moveBackward(self):
        self.camPos -= self.lookDir * self.moveSpeed
        self.updateLights()

    def moveLeft(self):
        self.camPos -= self.right * self.moveSpeed
        self.updateLights()

    def moveRight(self):
        self.camPos += self.right * self.moveSpeed
        self.updateLights()

    def rotateSystemXL(self):
        self.camRotGlobal[1] += self.globalRotSpeed

    def rotateSystemXR(self):
        self.camRotGlobal[1] -= self.globalRotSpeed

    def rotateSystemYL(self):
        self.camRotGlobal[2] += self.globalRotSpeed

    def rotateSystemYR(self):
        self.camRotGlobal[2] -= self.globalRotSpeed

    def rotateSystemZL(self):
        self.camRotGlobal[0] += self.globalRotSpeed

    def rotateSystemZR(self):
        self.camRotGlobal[0] -= self.globalRotSpeed

    def calcCameraDirections(self):
        self.lookDir = np.array([cos(self.camRot[1]) * sin(self.camRot[0]),
                                 sin(self.camRot[1]),
                                 cos(self.camRot[1]) * cos(self.camRot[0])])
        self.right = np.array(
            [sin(self.camRot[0] - pi / 2.0), 0, cos(self.camRot[0] - pi / 2.0)])

        self.up = np.cross(self.right, self.lookDir)

    def rotateCamera(self, mousePos, mousePosOld, mouseButtonState):
        dm = mousePos - mousePosOld

        if mouseButtonState[GLUT_LEFT_BUTTON] == GLUT_DOWN:
            self.camRotGlobal[1] += dm[0] * 0.1 * self.globalRotSpeed
            self.camRotGlobal[0] += dm[1] * 0.1 * self.globalRotSpeed
        elif mouseButtonState[GLUT_RIGHT_BUTTON] == GLUT_DOWN:
            self.camRot += dm * self.lookSpeed
            self.calcCameraDirections()

    def glLookAt(self):
        lookAt = self.camPos + self.lookDir
        gluLookAt(self.camPos[0], self.camPos[1], self.camPos[2],
                  lookAt[0], lookAt[1], lookAt[2],
                  self.up[0], self.up[1], self.up[2])

    def rotateSystem(self):
        glTranslatef(self.center[0], self.center[1], self.center[2])
        glRotatef(self.camRotGlobal[0], 1, 0, 0)
        glRotatef(self.camRotGlobal[1], 0, 1, 0)
        glRotatef(self.camRotGlobal[2], 0, 0, 1)
        glTranslatef(-self.center[0], -self.center[1], -self.center[2])
        self.updateLights()
