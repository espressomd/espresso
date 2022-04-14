.. _Running a simulation:

Running a simulation
====================

.. _Running es:

Running |es|
------------

Running a script
~~~~~~~~~~~~~~~~

|es| is implemented as a Python module. This means that you need to write a
python script for any task you want to perform with |es|. In this chapter,
the basic structure of the interface will be explained. For a practical
introduction, see the tutorials, which are also part of the
distribution. To use |es|, you need to import the espressomd module in your
Python script. To this end, the folder containing the python module
needs to be in the Python search path. The module is located in the
:file:`src/python` folder under the build directory. A convenient way to run
python with the correct path is to use the pypresso script located in
the build directory.

.. code-block:: bash

    ./pypresso simulation.py

The ``pypresso`` script is just a wrapper in order to expose the |es| python
module to the system's python interpreter by modifying the ``$PYTHONPATH``.
Please see the following chapter :ref:`Setting up the system` describing how
to actually write a simulation script for |es|.

Running an interactive notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running a Jupyter session requires using the ``ipypresso`` script, which
is also located in the build directory (its name comes from the IPython
interpreter, today known as Jupyter). To run the tutorials, you will need
to start a Jupyter session:

.. code-block:: bash

    cd doc/tutorials
    ../../ipypresso notebook

You may then browse through the different tutorial folders. Files whose name
ends with extension .ipynb can be opened in the browser. Click on the Run
button to execute the current block, or use the keyboard shortcut Shift+Enter.
If the current block is a code block, the ``In [ ]`` label to the left will
change to ``In [*]`` while the code is being executed, and become ``In [1]``
once the execution has completed. The number increments itself every time a
code cell is executed. This bookkeeping is extremely useful when modifying
previous code cells, as it shows which cells are out-of-date. It's also
possible to run all cells by clicking on the "Run" drop-down menu, then on
"Run All Below". This will change all labels to ``In [*]`` to show that the
first one is running, while the subsequent ones are awaiting execution.
You'll also see that many cells generate an output. When the output becomes
very long, Jupyter will automatically put it in a box with a vertical scrollbar.
The output may also contain static plots, dynamic plots and videos. It is also
possible to start a 3D visualizer in a new window, however closing the window
will exit the Python interpreter and Jupyter will notify you that the current
Python kernel stopped. If a cell takes too long to execute, you may interrupt
it with the stop button.

Solutions cells are created using the ``exercise2`` plugin from nbextensions.
To prevent solution code cells from running when clicking on "Run All", these
code cells need to be converted to Markdown cells and fenced with `````python``
and ```````.

To close the Jupyter session, go to the terminal where it was started and use
the keyboard shortcut Ctrl+C twice.

When starting a Jupyter session, you may see the following warning in the
terminal:

.. code-block:: none

    [TerminalIPythonApp] WARNING | Subcommand `ipython notebook` is deprecated and will be removed in future versions.
    [TerminalIPythonApp] WARNING | You likely want to use `jupyter notebook` in the future

This only means |es| was compiled with IPython instead of Jupyter. If Jupyter
is installed on your system, the notebook will automatically close IPython and
start Jupyter. To recompile |es| with Jupyter, provide ``cmake`` with the flag
``-DIPYTHON_EXECUTABLE=$(which jupyter)``.

You can find the official Jupyter documentation at
https://jupyter.readthedocs.io/en/latest/running.html

Running inside an IDE
~~~~~~~~~~~~~~~~~~~~~

You can use an integrated development environment (IDE) to develop and run |es|
scripts. Suitable IDEs are e.g. *Visual Studio Code* and *Spyder*. They can
provide a workflow superior to that of a standard text editor as they offer
useful features such as advanced code completion, debugging and analysis tools
etc. The following example shows how to setup |es| in *Visual Studio Code* on
Linux (tested with version 1.46.1). The process should be similar for every
Python IDE, namely the Python interpreter needs to be replaced.

The ``pypresso`` executable can be set as a custom Python interpreter inside VS
Code. |es| scripts can then be executed just like any other python script.
Inside VS Code, the Python extension needs to be installed. Next, click the
gear at the bottom left and choose *Settings*. Search for
``Default Interpreter Path`` and change the setting to the path to your
``pypresso`` executable, e.g.

.. code-block:: none

    ~/espresso/build/pypresso

After that, you can open scripts and execute them with the keyboard shortcut
Ctrl+F5.

Fig. :ref:`vs-code-figure` shows the VS Code interface with the interpreter
path set to ``pypresso``.

.. note:: You may need to set the path relative to your home directory, i.e. ``~/path/to/pypresso``.

.. _vs-code-figure:

.. figure:: figures/vs-code-settings.png
   :alt: Visual Studio Code interface with the default interpreter path set to the ``pypresso`` executable
   :width: 55.0%
   :align: center

   Visual Studio Code interface


.. _Debugging es:

Debugging |es|
--------------

Exceptional situations occur in every program.  If |es| crashes with a
segmentation fault, that means that there was a memory fault in the
simulation core which requires running the program in a debugger.  The
``pypresso`` executable file is actually not a program but a script
which sets the Python path appropriately and starts the Python
interpreter with your arguments.  Thus it is not possible to directly
run ``pypresso`` in a debugger.  However, we provide some useful
command line options for the most common tools.

.. code-block:: bash

     ./pypresso --tool <args>

where ``--tool`` can be any from the following table.  You can only
use one tool at a time.

+---------------------+----------------------------------------------+
| Tool                | Effect                                       |
+=====================+==============================================+
| ``--gdb``           | ``gdb --args python <args>``                 |
+---------------------+----------------------------------------------+
| ``--lldb``          | ``lldb -- python <args>``                    |
+---------------------+----------------------------------------------+
| ``--valgrind``      | ``valgrind --leak-check=full python <args>`` |
+---------------------+----------------------------------------------+
| ``--cuda-gdb``      | ``cuda-gdb --args python <args>``            |
+---------------------+----------------------------------------------+
| ``--cuda-memcheck`` | ``cuda-memcheck python <args>``              |
+---------------------+----------------------------------------------+


.. _CUDA:

CUDA
----

:py:meth:`~espressomd.cuda_init.CudaInitHandle()` command can be used to choose the GPU for all subsequent
GPU-computations. Note that due to driver limitations, the GPU cannot be
changed anymore after the first GPU-using command has been issued, for
example ``lbfluid``. If you do not choose the GPU manually before that,
CUDA internally chooses one, which is normally the most powerful GPU
available, but load-independent. ::

    system = espressomd.System(box_l=[1, 1, 1])
    dev = system.cuda_init_handle.device
    system.cuda_init_handle.device = dev

The first invocation in the sample above returns the id of the set graphics card, the second one sets the
device id.

.. _GPU Acceleration with CUDA:

GPU Acceleration with CUDA
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    Feature ``CUDA`` required

|es| is capable of GPU acceleration to speed up simulations.
Not every simulation method is parallelizable or profits from
GPU acceleration. Refer to :ref:`Available simulation methods`
to check whether your desired method can be used on the GPU.
In order to use GPU acceleration you need a NVIDIA GPU
and it needs to have at least compute capability 2.0.

For more information please check :class:`espressomd.cuda_init.CudaInitHandle`.

.. _List available CUDA devices:

List available CUDA devices
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to list available CUDA devices, you should call
:meth:`espressomd.cuda_init.CudaInitHandle.list_devices`::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> print(system.cuda_init_handle.list_devices())
    {0: 'GeForce RTX 2080', 1: 'GeForce GT 730'}

This method returns a dictionary containing
the device id as key and the device name as its value.

To get more details on the CUDA devices for each MPI node, call
:meth:`espressomd.cuda_init.CudaInitHandle.list_devices_properties`::

    >>> import pprint
    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> pprint.pprint(system.cuda_init_handle.list_devices_properties())
    {'seraue': {0: {'name': 'GeForce RTX 2080',
                    'compute_capability': (7, 5),
                    'cores': 46,
                    'total_memory': 8370061312},
                1: {'name': 'GeForce GT 730',
                    'compute_capability': (3, 5),
                    'cores': 2,
                    'total_memory': 1014104064}}}

.. _Selection of CUDA device:

Selection of CUDA device
~~~~~~~~~~~~~~~~~~~~~~~~

When you start ``pypresso`` your first GPU should be selected.
If you wanted to use the second GPU, this can be done
by setting :attr:`espressomd.cuda_init.CudaInitHandle.device` as follows::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> system.cuda_init_handle.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.
