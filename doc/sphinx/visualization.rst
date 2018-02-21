.. _Online-visualization:

Online-visualization
====================

With the python interface, |es| features two possibilities for
online-visualization:

#. Using the mlab module to drive *Mayavi, a "3D scientific data
   visualization and plotting in Python"*. Mayavi has a user-friendly
   GUI to specify the appearance of the output.
   Additional requirements:
   python module *mayavi*, VTK (package *python-vtk* for Debian/Ubuntu).
   Note that only VTK from version 7.0.0 and higher has Python 3
   support.

#. A direct rendering engine based on *pyopengl*. As it is developed for |es|, 
   it supports the visualization of several specific features like
   external forces or constraints. It has no GUI to setup the
   appearance, but can be adjusted by a large set of parameters.
   Additional requirements:
   python module *PyOpenGL*.

Both are not meant to produce high quality renderings, but rather to
debug your setup and equilibration process.

.. _General usage:

General usage
-------------

The recommended usage of both tools is similar: Create the visualizer of
your choice and pass it the ``espressomd.System()`` object. Then write
your integration loop in a separate function, which is started in a
non-blocking thread. Whenever needed, call ``update()`` to synchronize
the renderer with your system. Finally start the blocking visualization
window with ``start()``. See the following minimal code example::

    import espressomd 
    from espressomd import visualization 
    from threading import Thread

    system = espressomd.System() 
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    system.box_l = [10,10,10]

    system.part.add(pos = [1,1,1]) 
    system.part.add(pos = [9,9,9])

    #visualizer = visualization.mayaviLive(system) 
    visualizer = visualization.openGLLive(system)

    def main_thread(): 
        while True: 
            system.integrator.run(1)
            visualizer.update()

    t = Thread(target=main_thread) 
    t.daemon = True 
    t.start()
    visualizer.start()

.. _Common methods for openGL and mayavi:

Common methods for openGL and mayavi
------------------------------------

| :meth:`espressomd.visualization.mayaviLive.update()` 
| :meth:`espressomd.visualization.openGLLive.update()`

``update()`` synchronizes system and visualizer, handles keyboard events for
openGLLive.

| :meth:`espressomd.visualization.mayaviLive.start()` 
| :meth:`espressomd.visualization.openGLLive.start()`

``start()`` starts the blocking visualizer window. 
Should be called after a separate thread containing ``update()`` has been started.

| :meth:`espressomd.visualization.mayaviLive.register_callback()`
| :meth:`espressomd.visualization.openGLLive.register_callback()`

Registers the method ``callback()``, which is called every ``interval`` milliseconds. Useful for
live plotting (see sample script samples/python/visualization.py).

.. _Mayavi visualizer:

Mayavi visualizer
-----------------

The mayavi visualizer is created with the following syntax:

:class:`espressomd.visualization.mayaviLive()`

Required parameters:
    * `system`: The espressomd.System() object.
Optional keywords:
    * `particle_sizes`:
        * `"auto"` (default)`: The Lennard-Jones sigma value of the self-interaction is used for the particle diameter.
        * `callable`: A lambda function with one argument. Internally, the numerical particle type is passed to the lambda function to determine the particle radius.
        * `list`: A list of particle radii, indexed by the particle type.

.. _OpenGL visualizer:

OpenGL visualizer
-----------------

| :meth:`espressomd.visualization.openGLLive.run()` 

To visually debug your simulation, ``run()`` can be used to conveniently start 
an integration loop in a separate thread once the visualizer is initialized::

    import espressomd 
    from espressomd import visualization 

    system = espressomd.System() 
    system.cell_system.skin = 0.4
    system.time_step = 0.00001
    system.box_l = [10,10,10]

    system.part.add(pos = [1,1,1], v = [1,0,0]) 
    system.part.add(pos = [9,9,9], v = [0,1,0])

    visualizer = visualization.openGLLive(system, background_color = [1,1,1])
    visualizer.run(1)

:class:`espressomd.visualization.openGLLive()`

The optional keywords in ``**kwargs`` are used to adjust the appearance of the visualization.
The parameters have suitable default values for most simulations. 

Required parameters:
    * `system`: The espressomd.System() object.
Optional keywords:
    * `window_size`: Size of the visualizer window in pixels.
    * `name`: The name of the visualizer window.
    * `background_color`: RGB of the background.
    * `periodic_images`: Periodic repetitions on both sides of the box in xyzdirection.
    * `draw_box`: Draw wireframe boundaries.
    * `draw_axis`: Draws xyz system axes.
    * `quality_particles`: The number of subdivisions for particle spheres.
    * `quality_bonds`: The number of subdivisions for cylindrical bonds.
    * `quality_arrows`: The number of subdivisions for external force arrows.
    * `quality_constraints`: The number of subdivisions for primitive constraints.
    * `close_cut_distance`: The distance from the viewer to the near clipping plane.
    * `far_cut_distance`: The distance from the viewer to the far clipping plane.
    * `camera_position`: Initial camera position. `auto` (default) for shifted position in z-direction. 
    * `camera_target`: Initial camera target. `auto` (default) to look towards the system center.
    * `camera_right`: Camera right vector in system coordinates. Default is [1, 0, 0] 
    * `particle_sizes`:     
        * `auto` (default)`: The Lennard-Jones sigma value of the self-interaction is used for the particle diameter.
        * `callable`: A lambda function with one argument. Internally, the numerical particle type is passed to the lambda function to determine the particle radius.
        * `list`: A list of particle radii, indexed by the particle type.
    * `particle_coloring`:  
        * `auto` (default)`: Colors of charged particles are specified by particle_charge_colors, neutral particles by particle_type_colors
        * `charge`: Minimum and maximum charge of all particles is determined by the visualizer. All particles are colored by a linear interpolation of the two colors given by particle_charge_colors according to their charge.
        * `type`: Particle colors are specified by particle_type_colors, indexed by their numerical particle type.
    * `particle_type_colors`: Colors for particle types.
    * `particle_type_materials`: Materials of the particle types.
    * `particle_charge_colors`: Two colors for min/max charged particles.
    * `draw_constraints`: Enables constraint visualization. For simple constraints (planes, spheres and cylinders), OpenGL primitives are used. Otherwise, visualization by rasterization is used.
    * `rasterize_pointsize`: Point size for the rasterization dots.
    * `rasterize_resolution`: Accuracy of the rasterization.
    * `quality_constraints`: The number of subdivisions for primitive constraints.
    * `constraint_type_colors`: Colors of the constraints by type.
    * `constraint_type_materials`: Materials of the constraints by type.
    * `draw_bonds`: Enables bond visualization.
    * `bond_type_radius`: Radii of bonds by type.
    * `bond_type_colors`: Color of bonds by type.
    * `bond_type_materials`: Materials of bonds by type.
    * `ext_force_arrows`: Enables external force visualization.
    * `ext_force_arrows_scale`: Scale factor for external force arrows.
    * `drag_enabled`: Enables mouse-controlled particles dragging (Default`: False)
    * `drag_force`: Factor for particle dragging
    * `light_pos`: If `auto` (default) is used, the light is placed dynamically in the particle barycenter of the system. Otherwise, a fixed coordinate can be set.
    * `light_colors`: Three lists to specify ambient, diffuse and specular light colors.
    * `light_brightness`: Brightness (inverse constant attenuation) of the light. 
    * `light_size`: Size (inverse linear attenuation) of the light. If `auto` (default) is used, the light size will be set to a reasonable value according to the box size at start.
    * `spotlight_enabled`: If set to ``True`` (default), it enables a spotlight on the camera position pointing in look direction.
    * `spotlight_colors`: Three lists to specify ambient, diffuse and specular spotlight colors.
    * `spotlight_angle`: The spread angle of the spotlight in degrees (from 0 to 90).
    * `spotlight_brightness`: Brightness (inverse constant attenuation) of the spotlight. 
    * `spotlight_focus`: Focus (spot exponent) for the spotlight from 0 (uniform) to 128. 

.. _Colors and Materials:

Colors and Materials
~~~~~~~~~~~~~~~~~~~~

Colors for particles, bonds and constraints are specified by RGBA arrays.
Materials by an array for the ambient, diffuse, specular and shininess (ADSS)
components. To distinguish particle groups, arrays of RGBA or ADSS entries are
used, which are indexed circularly by the numerical particle type::

    # Particle type 0 is red, type 1 is blue (type 2 is red etc)..
    visualizer = visualization.openGLLive(system, 
                                          particle_coloring = 'type',
                                          particle_type_colors = [[1, 0, 0, 1],[0, 0, 1, 1]])

`particle_type_materials` lists the materials by type::

    # Particle type 0 is gold, type 1 is blue (type 2 is gold again etc).
    visualizer = visualization.openGLLive(system, 
                                          particle_coloring = 'type',
                                          particle_type_colors = [[1, 1, 1, 1],[0, 0, 1, 1]],
                                          particle_type_materials = [gold, bright])

Materials are stored in :attr:`espressomd.visualization.openGLLive().materials`. 

.. _Controls:

Controls
~~~~~~~~

The camera can be controlled via mouse and keyboard:

    * hold left button: rotate the system
    * hold right button: translate the system
    * hold middle button: zoom / roll
    * mouse wheel / key pair TG: zoom
    * WASD-Keyboard control (WS: move forwards/backwards, AD: move sidewards)
    * Key pairs QE, RF, ZC: rotate the system 

Additional input functionality for mouse and keyboard is possible by assigning
callbacks to specified keyboard or mouse buttons. This may be useful for
realtime adjustment of system parameters (temperature, interactions, particle
properties etc) of for demonstration purposes. The callbacks can be triggered
by a timer or keyboard input:: 

    def foo():
        print "foo"

    #Registers timed calls of foo()
    visualizer.register_callback(foo,interval=500)

    #Callbacks to control temperature 
    temperature = 1.0
    def increaseTemp():
            global temperature
            temperature += 0.1
            system.thermostat.set_langevin(kT=temperature, gamma=1.0)
            print "T =",system.thermostat.get_state()[0]['kT']

    def decreaseTemp():
        global temperature
        temperature -= 0.1

        if temperature > 0:
            system.thermostat.set_langevin(kT=temperature, gamma=1.0)
            print "T =",system.thermostat.get_state()[0]['kT']
        else:
            temperature = 0
            system.thermostat.turn_off()
            print "T = 0"

    #Registers input-based calls
    visualizer.keyboardManager.registerButton(KeyboardButtonEvent('t',KeyboardFireEvent.Hold,increaseTemp))
    visualizer.keyboardManager.registerButton(KeyboardButtonEvent('g',KeyboardFireEvent.Hold,decreaseTemp))

Further examples can be found in samples/python/billard.py or samples/python/visualization\_openGL.py.

.. _Dragging particles:

Dragging particles
~~~~~~~~~~~~~~~~~~

With the keyword ``drag_enabled`` set to ``True``, the mouse can be used to
exert a force on particles in drag direction (scaled by ``drag_force`` and the
distance of particle and mouse cursor). 

.. _Visualization example scripts:

Visualization example scripts
-----------------------------

Various example scripts can be found in the samples/python folder or in
some tutorials:

-  samples/python/visualization.py: LJ-Liquid with live plotting.

-  samples/python/visualization\_bonded.py: Sample for bond
   visualization.

-  samples/python/billard.py: Simple billard game including many
   features of the openGL visualizer.

-  samples/python/visualization\_openGL.py: Timer and keyboard callbacks
   for the openGL visualizer.

-  doc/tutorials/python/02-charged\_system/scripts/nacl\_units\_vis.py:
   Periodic NaCl crystal, see tutorial “Charged Systems”.

-  doc/tutorials/python/02-charged\_system/scripts/nacl\_units\_confined\_vis.py:
   Confined NaCl with interactively adjustable electric field, see
   tutorial “Charged Systems”.

-  doc/tutorials/python/08-visualization/scripts/visualization.py:
   LJ-Liquid visualization along with tutorial “Visualization”.

Finally, it is recommended to go through tutorial “Visualization” for
further code explanations. Also, the tutorial “Charged Systems” has two
visualization examples.
