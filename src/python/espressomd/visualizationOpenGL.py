from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
from math import *
import numpy as np
import os
import time

class openGLLive:

    def __init__(self, system, specs={}):

        #TODO: COMPLETE DESCRIPTION FOR SPECS

        #		'particle_coloring':		 'auto', 'charge', 'type'
        #		'particle_type_colors':		 [r_t1,g_t1,b_t1],[r_t2,g_t2,b_t2],... ]
        #		'particle_charge_colors':	 [[r_lowq,g_lowq,b_lowq],[r_highq,g_highq,b_highq]]
        #		'particle_sizes':		     'auto', 'type'
        #		'particle_type_sizes':			     [size_t1,size_t2,..]
        #		'ext_force_arrows': 		 True/False
        #		'window_size':				 [x,y]
        #		'background_color':			 [r,g,b]
        #		'update_fps':				 fps
        #	    'draw_bonds':				 True/False
        #       'bond_coloring':			 'type'
        #       'bond_type_radius':			 [r_t1,r_t2,...]
        #       'bond_type_colors':			 [r_t1,g_t1,b_t1],[r_t2,g_t2,b_t2],... ]
        #       'LB':						 True/False
        #		'light_pos':				 'auto', [x,y,z]
        #       'light_color':				 [r,g,b]
        #       'lightDecay':				 factor*[box_l] to light attenuation
        #       'particle_type_materials'

        #USER FRIENDLY DICT WITH VISUALIZATION SPECIFICATIONS

        self.specs = {
            'window_size':				  [800, 800],
            'name':						  'Espresso Visualization',
            'background_color':			  [0, 0, 0],
            'update_fps':				  30,
            'periodic_images':	   		  [0, 0, 0],
            'draw_box':		 	  		  True,
            'quality_spheres':            15,
            'quality_bonds':              15,
            'quality_arrows':             15,
            'close_cut_distance':         0.1,
            'far_cut_distance':           5,

            'particle_coloring':   		  'auto',
            'particle_sizes':			  'auto',
            'particle_type_colors':		  [[1, 1, 0, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],
            'particle_type_materials':	  [[0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1], [0.6, 1, 0.1]],
            'particle_charge_colors':	  [np.array([1, 0, 0, 1]), np.array([0, 1, 0, 1])],
            'particle_type_sizes':		  [1, 1, 1, 1, 1, 1, 1, ],

            'draw_bonds':				  True,
            'bond_coloring':			  'type',
            'bond_type_radius':		      [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
            'bond_type_colors':			  [[1, 1, 1, 1], [1, 0, 1, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 1, 0, 1], [1, 0.5, 0, 1], [0.5, 0, 1, 1]],

            'ext_force_arrows': 		  True,
            'ext_force_arrows_scale': 	  [1, 1, 1, 1, 1, 1, 1],

            'LB':						  False,

            'light_pos':				  'auto',
            'light_color':				  [0.8, 0.8, 0.8],
            'lightBrightness':   		  1.0,
            'lightSize':		   		  1.0,

            'dragEnabled':		   		  True,
            'dragForce':		   		  3.0
        }

        for prop in specs.iterkeys():
            if prop not in self.specs.iterkeys():
                raise ValueError(
                    prop + 'is not a valid visualization property')
            else:
                self.specs[prop] = specs[prop]

        self.system = system
        self.started = False
        self.keyboardManager = KeyboardManager()
        self.mouseManager = MouseManager()
        self.timers = []

    #CALLBACKS FOR THE MAIN THREAD
    def registerCallback(self, cb, interval=1000):
        self.timers.append((int(interval), cb))

    #THE BLOCKING START METHOD
    def start(self):
        self.initOpenGL()
        self.initEspressoVisualization()
        self.initCamera()
        self.initControls()
        self.initCallbacks()

        #POST DISPLAY WITH 60FPS
        def timed_update_redraw(data):
            glutPostRedisplay()
            glutTimerFunc(17, timed_update_redraw, -1)

        #PLACE LIGHT AT PARTICLE CENTER, DAMPED SPRING FOR SMOOTH POSITION CHANGE, CALL WITH 10FPS
        def timed_update_centerLight(data):
            ldt = 0.8
            cA = (self.particle_COM - self.smooth_light_pos) * \
                0.1 - self.smooth_light_posV * 1.8
            self.smooth_light_posV += ldt * cA
            self.smooth_light_pos += ldt * self.smooth_light_posV
            self.updateLightPos=True
            glutTimerFunc(100, timed_update_centerLight, -1)

        #AVERAGE PARTICLE COM ONLY EVERY 2sec
        def timed_update_particleCOM(data):
            self.particle_COM = np.average(self.particles['coords'], axis=0)
            glutTimerFunc(2000, timed_update_particleCOM, -1)

        self.started = True
        glutTimerFunc(17, timed_update_redraw, -1)
        if self.specs['light_pos'] == 'auto':
            glutTimerFunc(2000, timed_update_particleCOM, -1)
            glutTimerFunc(60, timed_update_centerLight, -1)
        #FOR MAC, BRING WINDOW TO FRONT
        os.system(
            '''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        #START THE BLOCKING MAIN LOOP
        glutMainLoop()

    #CALLED FROM ESPRESSO INTEGRATION LOOP
    #CHANGES OF ESPRESSO SYSTEM CAN ONLY HAPPEN HERE
    def update(self):
        if self.started:
            
            #IF CALLED TOO OFTEN, ONLY UPDATE WITH GIVEN FREQ
            self.elapsedTime += (time.time() - self.measureTimeBeforeIntegrate)
            if self.elapsedTime > 1.0 / self.specs['update_fps']:
                self.elapsedTime = 0
                self.updateParticles()
                #KEYBOARD CALLBACKS MAY CHANGE ESPRESSO SYSTEM PROPERTIES, ONLY SAVE TO CHANGE HERE
                self.keyboardManager.handleInput()

            self.measureTimeBeforeIntegrate = time.time()

            if self.triggerSetParticleDrag==True and self.dragId != -1:
                self.system.part[self.dragId].ext_force = self.dragExtForce
                self.triggerSetParticleDrag=False
            elif self.triggerResetParticleDrag==True and self.dragId != -1:
                self.system.part[self.dragId].ext_force = self.extForceOld
                self.triggerResetParticleDrag=False
                self.dragId = -1


    #GET THE PARTICLE DATA
    def updateParticles(self):
        self.particles = {'coords':  	self.system.part[:].pos,
                          'types':   	self.system.part[:].type,
                          'ext_forces': self.system.part[:].ext_force,
                          'charges': 	self.system.part[:].q
                          }

    #GET THE BOND DATA, SO FAR CALLED ONCE UPON INITIALIZATION
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

    #DRAW CALLED AUTOMATICALLY FROM GLUT DISPLAY FUNC
    def draw(self):

        if self.specs['LB']:
            self.drawLBVel()
        if self.specs['draw_box']:
            self.drawSystemBox()
        self.drawSystemParticles()

#		drawSphere(self.smooth_light_pos,0.5,[0,1.0,0,1.0],[1.,1.,1.])
#		drawSphere(self.particle_COM,0.5,[1.0,0,0,1.0],[1.,1.,1.])

        if self.specs['draw_bonds']:
            self.drawBonds()

    def drawSystemBox(self):
        box_l = self.system.box_l[0]
        b2 = box_l / 2.0
        drawBox([b2, b2, b2], box_l, [1 - self.specs['background_color'][0], 1 -
                                      self.specs['background_color'][1], 1 - self.specs['background_color'][2]])

    def drawSystemParticles(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        for pid in pIds:
            pos = coords[pid]
            q = self.particles['charges'][pid]
            ptype = self.particles['types'][pid]
            ext_f = self.particles['ext_forces'][pid]

            # Size: Lennard Jones Sigma,
            if self.specs['particle_sizes'] == 'auto':
                lj_sig = self.system.non_bonded_inter[ptype, ptype].lennard_jones.get_params()[
                    'sigma']
                radius = lj_sig / 2.0
                if radius == 0:
                    radius = self.sizeByType(ptype)

            elif self.specs['particle_sizes'] == 'type':
                radius = self.sizeByType(ptype)

            material = self.materialByType(ptype)

            if self.specs['particle_coloring'] == 'id':
                color = self.IdToColorf(pid)
                glColor(color)
            elif self.specs['particle_coloring'] == 'auto':
                # Color auto: Charge then Type
                if q != 0:
                    color = self.colorByCharge(q)
                else:
                    color = self.colorByType(ptype)
            elif self.specs['particle_coloring'] == 'charge':
                color = self.colorByCharge(q)
            elif self.specs['particle_coloring'] == 'type':
                color = self.colorByType(q)

            drawSphere(pos, radius, color, material, self.specs['quality_spheres'])
            for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                    for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                        if imx != 0 or imy != 0 or imz != 0:
                            redrawSphere(pos + (imx * self.imPos[0]+imy*self.imPos[1]+imz*self.imPos[2]), radius, self.specs['quality_spheres'])

            if self.specs['ext_force_arrows'] or pid == self.dragId:
                if ext_f[0] != 0 or ext_f[1] != 0 or ext_f[2] != 0:
                    if pid == self.dragId:
                        sc = 1
                    else:
                        sc = self.extForceArrowScaleByType(ptype)
                    if sc > 0:
                        drawArrow(pos, np.array(ext_f) * sc, 0.25, [1, 1, 1], self.specs['quality_arrows'])

    def drawBonds(self):
        coords = self.particles['coords']
        pIds = range(len(coords))
        b2 = self.system.box_l[0] / 2.0
        box_l2_sqr = pow(b2, 2.0)
        for b in self.bonds:
            if self.specs['bond_coloring'] == 'type':
                col = self.bondColorByType(b[2])
            radius = self.bondRadiusByType(b[2])
            d = coords[b[0]] - coords[b[1]]
            bondLen_sqr = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]

            if bondLen_sqr < box_l2_sqr:
                drawCylinder(coords[b[0]], coords[b[1]], radius, col, self.specs['quality_bonds'])
                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                drawCylinder(coords[b[0]] + im * self.imPos[dim], coords[b[1]] + im * self.imPos[dim], radius, col, self.specs['quality_bonds'])
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
                drawCylinder(coords[b[0]], s0, radius, col, self.specs['quality_bonds'])
                drawCylinder(coords[b[1]], s1, radius, col, self.specs['quality_bonds'])

                for imx in range(-self.specs['periodic_images'][0], self.specs['periodic_images'][0] + 1):
                    for imy in range(-self.specs['periodic_images'][1], self.specs['periodic_images'][1] + 1):
                        for imz in range(-self.specs['periodic_images'][2], self.specs['periodic_images'][2] + 1):
                            if imx != 0 or imy != 0 or imz != 0:
                                drawCylinder(coords[b[0]] + im * self.imPos[dim], s0 + im * self.imPos[dim], radius, col, self.specs['quality_bonds'])
                                drawCylinder(coords[b[1]] + im * self.imPos[dim], s1 + im * self.imPos[dim], radius, col, self.specs['quality_bonds'])

    #HELPER TO DRAW PERIODIC BONDS
    def isInsideBox(self, p):
        eps = 1e-5
        for i in range(3):
            if p[i] < -eps or p[i] > eps + self.system.box_l[i]:
                return False
        return True

    #VOXELS FOR LB VELOCITIES
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

    #USE MODULO IF THERE ARE MORE PARTICLE TYPES THAN TYPE DEFINITIONS FOR COLORS, MATERIALS ETC..
    def extForceArrowScaleByType(self, btype):
        return self.specs['ext_force_arrows_scale'][btype % len(self.specs['ext_force_arrows_scale'])]

    def materialByType(self, btype):
        return self.specs['particle_type_materials'][btype % len(self.specs['particle_type_materials'])]

    def bondColorByType(self, btype):
        return self.specs['bond_type_colors'][btype % len(self.specs['bond_type_colors'])]

    def bondRadiusByType(self, btype):
        return self.specs['bond_type_radius'][btype % len(self.specs['bond_type_radius'])]

    def sizeByType(self, ptype):
        return self.specs['particle_type_sizes'][ptype % len(self.specs['particle_type_sizes'])]

    def colorByType(self, ptype):
        return self.specs['particle_type_colors'][ptype % len(self.specs['particle_type_colors'])]

    #FADE PARTICE CHARGE COLOR FROM WHITE (q=0) to PLUSCOLOR (q=q_max) RESP MINUSCOLOR (q=q_min)
    def colorByCharge(self, q):
        if q < 0:
            c = 1.0 * q / self.minq
            return self.specs['particle_charge_colors'][0] * c + (1 - c) * np.array([1, 1, 1, 1])
        else:
            c = 1.0 * q / self.maxq

            return self.specs['particle_charge_colors'][1] * c + (1 - c) * np.array([1, 1, 1, 1])

    #ON INITIALIZATION, CHECK q_max/q_min
    def updateChargeColorRange(self):
        if len(self.particles['charges'][:])>0:
            self.minq = min(self.particles['charges'][:])
            self.maxq = max(self.particles['charges'][:])

    # INITS FOR GLUT FUNCTIONS
    def initCallbacks(self):
        # OpenGl Callbacks
        def display():
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            glLoadIdentity()

            self.camera.glLookAt()
            self.camera.rotateSystem()

            if self.updateLightPos:
                self.setLightPos()
                self.updateLightPos=False

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

        #CALLED ION WINDOW POSITION/SIZE CHANGE
        def reshapeWindow(w, h):
            glViewport(0, 0, w, h)
            glMatrixMode(GL_PROJECTION)
            glLoadIdentity()
            box_diag = pow(pow(self.system.box_l[0], 2) + pow(self.system.box_l[1], 2) + pow(self.system.box_l[1], 2), 0.5)
            gluPerspective(40, 1.0*w/h, self.specs['close_cut_distance'], self.specs['far_cut_distance'] * box_diag)
            glMatrixMode(GL_MODELVIEW)
            self.specs['window_size'][0] =1.0*w
            self.specs['window_size'][1] =1.0*h
            #glPushMatrix()

        #TIMERS FOR registerCallback
        def dummyTimer(index):
            self.timers[index][1]()
            glutTimerFunc(self.timers[index][0], dummyTimer, index)

        glutDisplayFunc(display)
        glutMouseFunc(mouse)
        glutKeyboardFunc(keyboardDown)
        glutKeyboardUpFunc(keyboardUp)
        glutReshapeFunc(reshapeWindow)
        #TODO: ZOOM WITH MOUSEWHEEL
        # glutMouseWheelFunc(mouseWheel);
        glutMotionFunc(motion)

        index=0
        for t in self.timers:
            glutTimerFunc(t[0], dummyTimer, index)
            index+=1

    #CLICKED ON PARTICLE: DRAG; CLICKED ON BACKGROUND: CAMERA
    def mouseMotion(self, mousePos, mousePosOld):

        if self.dragId != -1:
            ppos = self.particles['coords'][self.dragId]
            viewport = glGetIntegerv(GL_VIEWPORT)
            mouseWorld = gluUnProject(mousePos[0], viewport[3] - mousePos[1], self.depth)

            self.dragExtForce = self.specs['dragForce'] * (np.asarray(mouseWorld) - np.array(ppos))
            self.triggerSetParticleDrag = True
            #self.system.part[self.dragId].ext_force = f
        else:
            self.camera.rotateCamera(mousePos, mousePosOld)

    #DRAW SCENE AGAIN WITHOUT LIGHT TO IDENTIFY PARTICLE ID BY PIXEL COLOR
    def setParticleDrag(self, pos, pos_old):

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        self.camera.glLookAt()
        self.camera.rotateSystem()

        oldColMode = self.specs['particle_coloring']
        self.specs['particle_coloring'] = 'id'
        glDisable(GL_LIGHTING)
        self.drawSystemParticles()
        viewport = glGetIntegerv(GL_VIEWPORT)

        readPixel = glReadPixelsui(
            pos[0], viewport[3] - pos[1], 1, 1, GL_RGB, GL_FLOAT)[0][0]
        depth = glReadPixelsf(
            pos[0], viewport[3] - pos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)[0][0]
        pid = self.fcolorToId(readPixel)
        print "Selected Particle ", pid

        self.dragId = pid
        if pid != -1:
            self.dragPosInitial = self.particles['coords'][self.dragId]
            self.extForceOld = self.particles['ext_forces'][self.dragId][:]
            self.depth = depth
        self.specs['particle_coloring'] = oldColMode
        glEnable(GL_LIGHTING)

    def resetParticleDrag(self, pos, pos_old):
        if self.dragId != -1:
            self.triggerResetParticleDrag = True

    def IdToColorf(self, pid):
        pid += 1
        return [int(pid / (256 * 256)) / 255.0, int((pid % (256 * 256)) / 256) / 255.0, (pid % 256) / 255.0, 1.0]

    def fcolorToId(self, fcol):
        return 256 * 256 * int(fcol[0] * 255) + 256 * int(fcol[1] * 255) + int(fcol[2] * 255) - 1

    #ALL THE INITS
    def initEspressoVisualization(self):
        self.maxq = 0
        self.minq = 0

        self.dragId = -1
        self.dragPosInitial = []
        self.extForceOld = []
        self.dragExtForceOld = []
        self.triggerResetParticleDrag=False
        self.triggerSetParticleDrag=False

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
        self.updateParticles()
        self.updateChargeColorRange()
        self.updateBonds()

    #BOX PLANES (NORMAL, ORIGIN) FOR PERIODIC BONDS
    def boxSizeDependence(self):
        self.box_n = [np.array([1, 0, 0]), np.array([0, 1, 0]), np.array(
            [0, 0, 1]), np.array([-1, 0, 0]), np.array([0, -1, 0]), np.array([0, 0, -1])]
        self.box_p = [np.array([0, 0, 0]), np.array([0, 0, 0]), np.array([0, 0, 0]), np.array(
            self.system.box_l), np.array(self.system.box_l), np.array(self.system.box_l)]

    #DEFAULT CONTROLS
    def initControls(self):
        self.mouseManager.registerButton(MouseButtonEvent(
            None, MouseFireEvent.FreeMotion, self.mouseMotion))
        if self.specs['dragEnabled']:
            self.mouseManager.registerButton(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonPressed, self.setParticleDrag))
            self.mouseManager.registerButton(MouseButtonEvent(
                GLUT_LEFT_BUTTON, MouseFireEvent.ButtonReleased, self.resetParticleDrag))

        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'w', KeyboardFireEvent.Hold, self.camera.moveForward))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            's', KeyboardFireEvent.Hold, self.camera.moveBackward))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'a', KeyboardFireEvent.Hold, self.camera.moveLeft))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'd', KeyboardFireEvent.Hold, self.camera.moveRight))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'e', KeyboardFireEvent.Hold, self.camera.rotateSystemXR))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'q', KeyboardFireEvent.Hold, self.camera.rotateSystemXL))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'z', KeyboardFireEvent.Hold, self.camera.rotateSystemYR))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'c', KeyboardFireEvent.Hold, self.camera.rotateSystemYL))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'r', KeyboardFireEvent.Hold, self.camera.rotateSystemZR))
        self.keyboardManager.registerButton(KeyboardButtonEvent(
            'f', KeyboardFireEvent.Hold, self.camera.rotateSystemZL))
    
    #ASYNCHRONOUS PARALLEL CALLS OF glLight CAUSES SEG FAULTS, SO ONLY CHANGE LIGHT AT CENTRAL display METHOD AND TRIGGER CHANGES
    def setLightPos(self): 
        if self.specs['light_pos'] == 'auto':
            glLightfv(GL_LIGHT0, GL_POSITION, [self.smooth_light_pos[0], self.smooth_light_pos[1], self.smooth_light_pos[2], 0.6])
        else:
            glLightfv(GL_LIGHT0, GL_POSITION, self.specs['light_pos'])

    def triggerLightPosUpdate(self):
        self.updateLightPos=True

    def initCamera(self):
        bl = self.system.box_l[0]
        bl2 = bl / 2.0
        box_center = np.array([bl2, bl2, bl2])
        self.camera = Camera(camPos=np.array([bl * 1.3, bl * 1.3, bl * 2.5]), camRot=np.array([3.55, -0.4]), center=box_center, updateLights=self.triggerLightPosUpdate)
        self.smooth_light_pos = np.copy(box_center)
        self.smooth_light_posV = np.array([0.0, 0.0, 0.0])
        self.particle_COM = np.copy(box_center)
        self.updateLightPos=True

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

        glLineWidth(2.0)
        glutIgnoreKeyRepeat(1)
        
        # setup lighting
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHTING)
        if self.specs['light_pos'] != 'auto':
            glLightfv(GL_LIGHT0, GL_POSITION, self.specs['light_pos'])
        else:
            glLightfv(GL_LIGHT0, GL_POSITION, self.system.box_l * 0.5)

        glLightfv(GL_LIGHT0, GL_DIFFUSE, self.specs['light_color'])

        glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION,
                 0.7 / self.specs['lightBrightness'])
        glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, self.system.box_l[
                 0] / 67. * 0.005 / self.specs['lightSize'])
        glEnable(GL_LIGHT0)

#END OF MAIN CLASS 

#OPENGL DRAW WRAPPERS
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

def drawBox(pos, size, color):
    setSolidMaterial(color[0], color[1], color[2], 1, 2, 1)
    glPushMatrix()
    glTranslatef(pos[0], pos[1], pos[2])
    glutWireCube(size)
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

def drawCylinder(posA, posB, radius, color, quality):
    setSolidMaterial(color[0], color[1], color[2])
    glPushMatrix()

    d = posB - posA
    v = np.linalg.norm(d)
    ax = 57.2957795 * acos(d[2] / v)
    if d[2] < 0.0:
        ax = -ax
    rx = -d[1] * d[2]
    ry = d[0] * d[2]

    #angle,t,length = calcAngle(d)
    length = np.linalg.norm(d)
    glTranslatef(posA[0], posA[1], posA[2])
    glRotatef(ax, rx, ry, 0.0)
    quadric = gluNewQuadric()
    gluCylinder(quadric, radius, radius, length, quality, quality)
    # glutSolidCylinder(radius,radius,length,10,10)
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

    #angle,t,length = calcAngle(d)
    glTranslatef(pos[0], pos[1], pos[2])
    glRotatef(ax, rx, ry, 0.0)
    # glRotatef(angle,t[0],t[1],t[2]);
    quadric = gluNewQuadric()
    gluCylinder(quadric, radius, radius, v, quality, quality)

    glPopMatrix()
    e = pos + d
    glTranslatef(e[0], e[1], e[2])
    # glRotatef(angle,t[0],t[1],t[2]);
    glRotatef(ax, rx, ry, 0.0)

    glutSolidCone(radius * 3, 3, quality, quality)
    glPopMatrix()


#MOUSE EVENT MANAGER
class MouseFireEvent:

    ButtonPressed = 0
    FreeMotion = 1
    ButtonMotion = 2
    ButtonReleased = 3


class MouseButtonEvent:

    def __init__(self, button, fireEvent, callback):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback


class MouseManager:

    def __init__(self):
        self.mousePos = np.array([0, 0])
        self.mousePosOld = np.array([0, 0])
        self.mouseEventsPressed = []
        self.mouseEventsFreeMotion = []
        self.mouseEventsButtonMotion = []
        self.mouseEventsReleased = []

    def registerButton(self, mouseEvent):
        if mouseEvent.fireEvent == MouseFireEvent.ButtonPressed:
            self.mouseEventsPressed.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.FreeMotion:
            self.mouseEventsFreeMotion.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonMotion:
            self.mouseEventsButtonMotion.append(mouseEvent)
        elif mouseEvent.fireEvent == MouseFireEvent.ButtonReleased:
            self.mouseEventsReleased.append(mouseEvent)

    def mouseClick(self, button, state, x, y):
        self.mousePosOld = self.mousePos
        self.mousePos = np.array([x, y])

        for me in self.mouseEventsPressed:
            if me.button == button and state == GLUT_DOWN:
                me.callback(self.mousePos, self.mousePosOld)
        for me in self.mouseEventsReleased:
            if me.button == button and state == GLUT_UP:
                me.callback(self.mousePos, self.mousePosOld)

    def mouseMove(self, x, y):
        self.mousePosOld = self.mousePos
        self.mousePos = np.array([x, y])

        for me in self.mouseEventsFreeMotion:
            me.callback(self.mousePos, self.mousePosOld)
#		for me in self.mouseEventsButtonMotion:
#			if me.button == button:
#				me.callback(self.mousePos,self.mousePosOld)


#KEYBOARD EVENT MANAGER
class KeyboardFireEvent:

    Pressed = 0
    Hold = 1
    Released = 2


class KeyboardButtonEvent:

    def __init__(self, button, fireEvent, callback):
        self.button = button
        self.fireEvent = fireEvent
        self.callback = callback


class KeyboardManager:

    def __init__(self):
        self.pressedKeys = set([])
        self.keyStateOld = {}
        self.keyState = {}
        self.buttonEventsPressed = []
        self.buttonEventsHold = []
        self.buttonEventsReleased = []

    def registerButton(self, buttonEvent):
        if buttonEvent.fireEvent == KeyboardFireEvent.Pressed:
            self.buttonEventsPressed.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Hold:
            self.buttonEventsHold.append(buttonEvent)
        elif buttonEvent.fireEvent == KeyboardFireEvent.Released:
            self.buttonEventsReleased.append(buttonEvent)

    def handleInput(self):

        removeKeys = set([])
        for b in self.pressedKeys:
            if self.keyStateOld[b] == 0 and self.keyState[b] == 1:
                for be in self.buttonEventsPressed:
                    if be.button == b:
                        be.callback()
                for be in self.buttonEventsHold:
                    if be.button == b:
                        be.callback()
#				print 'Key',b,'Pressed'

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 1:
                for be in self.buttonEventsHold:
                    if be.button == b:
                        be.callback()
#				print 'Key',b,'Hold'

            elif self.keyStateOld[b] == 1 and self.keyState[b] == 0:
                for be in self.buttonEventsReleased:
                    if be.button == b:
                        be.callback()
#				print 'Key',b,'Released'
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

#CAMERA
class Camera:

    def __init__(self, camPos=np.array([0, 0, 1]), camRot=np.array([pi, 0]), moveSpeed=3, rotSpeed=0.001, globalRotSpeed=3, center=np.array([0, 0, 0]), updateLights=None):
        self.moveSpeed = moveSpeed
        self.lookSpeed = rotSpeed
        self.globalRotSpeed = globalRotSpeed
        self.camPos = camPos
        self.camRot = camRot
        self.center = center
        self.camRotGlobal = np.array([0, 0, 0])
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

    def rotateCamera(self, mousePos, mousePosOld):
        self.camRot += (mousePos - mousePosOld) * self.lookSpeed
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
