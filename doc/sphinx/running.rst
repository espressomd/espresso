.. _Running a simulation:

Running a simulation
====================

|es| is implemented as a Python module. This means that you need to write a
python script for any task you want to perform with |es|. In this chapter,
the basic structure of the interface will be explained. For a practical
introduction, see the tutorials, which are also part of the distribution.

Most users should consider building and then installing |es| locally.
In this way, |es| behaves like any regular Python package and will
be recognized by the Python interpreter and Jupyter notebooks.

Most developers prefer the ``pypresso`` resp. ``ipypresso`` wrapper scripts,
which export the build folder into the ``$PYTHONPATH`` environment variable
and then call ``python`` resp. ``jupyter``. They also introduce extra command
line options to help developers run simulations inside a debugger.
Command line examples in this chapter use the wrapper scripts instead of the
Python and Jupyter programs, although they are perfectly interchangeable
when not using a debugger.

.. _Running es:

Running |es|
------------

Running a script
~~~~~~~~~~~~~~~~

To use |es|, you need to import the ``espressomd`` module in your
Python script. To this end, the folder containing the python module
needs to be in the Python search path. The module is located in the
:file:`src/python` folder under the build directory.

A convenient way to run Python with the correct path is to use the
``pypresso`` script located in the build directory:

.. code-block:: bash

    ./pypresso simulation.py

The ``pypresso`` script is just a wrapper in order to expose the |es| python
module to the system's Python interpreter by modifying the ``$PYTHONPATH``.
If you have installed |es| from a Linux package manager that doesn't provide
the ``pypresso`` script, you will need to modify the ``$PYTHONPATH`` and
possibly the ``$LD_LIBRARY_PATH`` too, depending on which symbols are missing.

The next chapter, :ref:`Setting up the system`, will explain in more details
how to write a simulation script for |es|. If you don't have any script,
simply call one of the files listed in section :ref:`Sample scripts`.

Using the console
~~~~~~~~~~~~~~~~~

Since |es| can be manipulated like any other Python module, it is possible
to interact with it in a Python interpreter. Simply run the ``pypresso``
script without arguments to start a Python session:

.. code-block:: bash

    ./pypresso

Likewise, a Jupyter console can be started with the ``ipypresso`` script,
which is also located in the build directory:

.. code-block:: bash

    ./ipypresso console

The name comes from the IPython interpreter, today known as Jupyter.

Interactive notebooks
~~~~~~~~~~~~~~~~~~~~~

Tutorials are available as notebooks, i.e. they consist of a ``.ipynb``
file which contains both the source code and the corresponding explanations.
They can be viewed, changed and run interactively. To generate the tutorials
in the build folder, do:

.. code-block:: bash

    make tutorials

The tutorials contain solutions hidden with the ``exercise2`` NB extension.
Since this extension is only available for Jupyter Notebook, JupyterLab
users need to convert the tutorials:

.. code-block:: bash

    for f in doc/tutorials/*/*.ipynb; do
      ./pypresso doc/tutorials/convert.py exercise2 --to-jupyterlab ${f}
    done

Likewise, VS Code Jupyter users need to convert the tutorials:

.. code-block:: bash

    for f in doc/tutorials/*/*.ipynb; do
      ./pypresso doc/tutorials/convert.py exercise2 --to-vscode-jupyter ${f}
    done

To interact with notebooks, move to the directory containing the tutorials
and call the ``ipypresso`` script to start a local Jupyter session.

For Jupyter Notebook and IPython users:

.. code-block:: bash

    cd doc/tutorials
    ../../ipypresso notebook

For JupyterLab users:

.. code-block:: bash

    cd doc/tutorials
    ../../ipypresso lab

For VS Code Jupyter users, no action is needed if ``pypresso`` was set as
the interpreter path (see details in :ref:`Running inside an IDE`).

You may then browse through the different tutorial folders. Files whose name
ends with extension ``.ipynb`` can be opened in the browser. Click on the Run
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
``-D IPYTHON_EXECUTABLE=$(which jupyter)``.

You can find the official Jupyter documentation at
https://jupyter.readthedocs.io/en/latest/running.html

.. _Running inside an IDE:

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

.. _Running in the cloud:

Running in the cloud
~~~~~~~~~~~~~~~~~~~~

A `Gitpod <https://gitpod.io>`__ config file is provided to automatically
build |es| in its default configuration (`direct link
<https://gitpod.io/#https://github.com/espressomd/espresso>`__), which is
sufficient to run most tutorials. The Gitpod workspace can be accessed from
the `terminal via SSH <https://www.gitpod.io/docs/configure/ssh>`__ or from
a `web browser <https://www.gitpod.io/docs/configure/browser-settings>`__,
which uses the VS Code IDE.

To execute the tutorials, choose a Jupyter backend:

* VS Code Jupyter: navigate to ``ESPRESSO/build/doc/tutorials`` in the
  project tree and open the notebook files; if the kernel drop-down menu
  doesn't offer ``build/pypresso`` as a kernel, restart the VS Code IDE:
  quit the workspace by closing the browser tab, re-open the tab and
  click ``espressomd-espresso-...`` in the popup to restart the IDE
  (don't click on the green button "New Workspace")

* Jupyter Notebook:

  .. code-block:: bash

      cd ${GITPOD_REPO_ROOT}/build/doc/tutorials
      ../../ipypresso notebook --NotebookApp.allow_origin="$(gp url 8888)" \
          --port=8888 --no-browser

* JupyterLab:

  .. code-block:: bash

      cd ${GITPOD_REPO_ROOT}/build/doc/tutorials
      ../../ipypresso lab --NotebookApp.allow_origin="$(gp url 8888)" \
          --port=8888 --no-browser

For both Jupyter Notebook and JupyterLab, a notification will appear and say
that a new port 8888 has been made available. Click the orange "Make public"
button to open that port and then Ctrl+click one of the urls in the terminal
output to open the Jupyter backed in a pop-up window.

To start a workspace from a specific branch, use a link in the following form:
``https://gitpod.io/#https://github.com/user_name/espresso/tree/branch_name``,
where ``user_name`` and ``branch_name`` need to be adapted.


.. _Parallel computing:

Parallel computing
------------------

Many algorithms in |es| are designed to work with multiple MPI ranks.
However, not all algorithms benefit from MPI parallelization equally.
Several algorithms only use MPI rank 0 (e.g. :ref:`Reaction methods`).
|es| should work with most MPI implementations on the market;
see the :term:`MPI installation requirements <MPI>` for details.

.. _General syntax:

General syntax
~~~~~~~~~~~~~~

To run a simulation on several MPI ranks, for example 4, simply invoke
the ``pypresso`` script with the following syntax:

.. code-block:: bash

    mpiexec -n 4 ./pypresso simulation.py

The cell system is automatically split among the MPI ranks, and data
is automatically gathered on the main rank, which means a regular |es|
script can be executed in an MPI environment out-of-the-box. The number
of MPI ranks can be accessed via the system ``n_nodes`` state property.
The simulation box partition is controlled by the cell system
:attr:`~espressomd.cell_system.CellSystem.node_grid` property.
By default, MPI ranks are assigned in decreasing order, e.g. on 6 MPI ranks
``node_grid`` is ``[3, 2, 1]``. It is possible to re-assign the ranks by
changing the value of the ``node_grid`` property, however a few algorithms
(such as FFT-based electrostatic methods) only work for the default
partitioning scheme where values must be arranged in decreasing order.

::

    # get the number of ranks
    print(system.cell_system.get_state()["n_nodes"])
    # re-assign the ranks
    system.cell_system.node_grid = [2, 1, 3]
    system.cell_system.node_grid = [6, 1, 1]

There are alternative ways to invoke MPI on ``pypresso``, but they share
similar options. The number after the ``-n`` option is the number of ranks,
which needs to be inferior or equal to the number of *physical* cores on the
workstation. Command ``nproc`` displays the number of *logical* cores on the
workstation. For architectures that support hyperthreading, the number of
logical cores is an integer multiple of the number of physical cores,
usually 2. Therefore on a hyperthreaded workstation with 32 cores,
at most 16 cores can be used without major performance loss, unless
extra arguments are passed to the ``mpiexec`` program.

On cluster computers, it might be necessary to load the MPI library with
``module load openmpi`` or similar.

.. _Performance gain:

Performance gain
~~~~~~~~~~~~~~~~

Simulations executed in parallel with run faster, however the runtime
won't decrease linearly with the number of MPI ranks. MPI-parallel
simulations introduce several sources of overhead and latency:

* overhead of serializing, communicating and deserializing data structures
* extra calculations in the LB halo
* extra calculations in the ghost shell
  (see section :ref:`Internal particle organization` for more details)
* latency due to blocking communication (i.e. a node remains idle
  while waiting for a message from another node)
* latency due to blocking data collection for GPU
  (only relevant for GPU methods)
* latency due to context switching
* latency due to memory bandwidth

While good performance can be achieved up to 32 MPI ranks, allocating more
than 32 ranks to a simulation will not always lead to significantly improved
run times. The performance gain is highly sensitive to the algorithms used
by the simulation, for example GPU methods rarely benefit from more than
8 MPI ranks. Performance is also affected by the number of features enabled
at compile time, even when these features are not used by the simulation;
do not hesitate to remove all features not required by the
simulation script and rebuild |es| for optimal performance.

Benchmarking is often the best way to determine the optimal number of MPI
ranks for a given simulation setup. Please refer to the wiki chapter on
`benchmarking <https://github.com/espressomd/espresso/wiki/Development#Benchmarking>`__
for more details.

Runtime speed-up is not the only appeal of MPI parallelization. Another
benefit is the possibility to distribute a calculation over multiple
compute nodes in clusters and high-performance environments, and therefore
split the data structures over multiple machines. This becomes necessary
when running simulations with millions of particles, as the memory
available on a single compute node would otherwise saturate.

.. _Communication model:

Communication model
~~~~~~~~~~~~~~~~~~~

|es| was originally designed for the "flat" model of communication:
each MPI rank binds to a logical CPU core. This communication model
doesn't fully leverage shared memory on recent CPUs, such as `NUMA
architectures <https://en.wikipedia.org/wiki/Non-uniform_memory_access>`__,
and |es| currently doesn't support the hybrid
MPI+\ `OpenMP <https://www.openmp.org>`__ programming model.

The MPI+CUDA programming model is supported, although only one GPU can be
used for the entire simulation. As a result, a blocking *gather* operation
is carried out to collect data from all ranks to the main rank, and a
blocking *scatter* operation is carried out to transfer the result of the
GPU calculation from the main rank back to all ranks. This latency limits
GPU-acceleration to simulations running on fewer than 8 MPI ranks.
For more details, see section :ref:`GPU acceleration`.

.. _The MPI callbacks framework:

The MPI callbacks framework
"""""""""""""""""""""""""""

When starting a simulation with :math:`n` MPI ranks, |es| will internally
use MPI rank :math:`0` as the head node (also referred to as the "main rank")
and MPI ranks :math:`1` to :math:`n-1` as worker nodes. The Python interface
interacts only with the head node, and the head node forwards the information
to the worker nodes.

To put it another way, all worker nodes are idle until the user calls
a function that is designed to run in parallel,
in which case the head node calls the corresponding core function
and sends a request on the worker nodes to call the same core function.
The request can be a simple collective call, or a collective call with a
reduction if the function returns a value. The reduction can either:

- combine the :math:`n` results via a mathematical operation
  (usually a summation or a multiplication)
- discard the result of the :math:`n-1` worker nodes; this is done when
  all ranks return the same value, or when the calculation can only be
  carried out on the main rank but requires data from the other ranks
- return the result of one rank when the calculation can only be carried out
  by a specific rank; this is achieved by returning an *optional*, which
  contains a value on the rank that has access to the information necessary
  to carry out the calculation, while the other :math:`n-1` ranks return
  an empty optional

For more details on this framework, please refer to the Doxygen documentation
of the the C++ core file :file:`MpiCallbacks.hpp`.


.. _GPU acceleration:

GPU acceleration
----------------

.. _CUDA acceleration:

CUDA acceleration
~~~~~~~~~~~~~~~~~

.. note::
    Feature ``CUDA`` required

|es| is capable of delegating work to the GPU to speed up simulations.
Not every simulation method profits from GPU acceleration.
Refer to :ref:`Available simulation methods`
to check whether your desired method can be used on the GPU.
In order to use GPU acceleration you need a NVIDIA GPU
and it needs to have at least compute capability 2.0.
For more details, please refer to the installation section
:ref:`Nvidia GPU acceleration`.

For more information please check :class:`espressomd.cuda_init.CudaInitHandle`.

.. _List available devices:

List available devices
""""""""""""""""""""""

To list available CUDA devices, call
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

.. _Select a device:

Select a device
"""""""""""""""

When you start ``pypresso``, the first GPU should be selected.
If you wanted to use the second GPU, this can be done
by setting :attr:`espressomd.cuda_init.CudaInitHandle.device` as follows::

    >>> import espressomd
    >>> system = espressomd.System(box_l=[1, 1, 1])
    >>> system.cuda_init_handle.device = 1

Setting a device id outside the valid range or a device
which does not meet the minimum requirements will raise
an exception.


.. _Instrumentation:

Instrumentation
---------------

.. _Debugging:

Debugging
~~~~~~~~~

Exceptional situations occur in every program. If |es| crashes with a
fatal error, it is necessary to use a debugger to investigate the issue.
The tool should be chosen depending on the nature of the bug.
Most fatal errors fall into one of these categories:

* segmentation fault: typically due to uninitialized pointers, dangling
  pointers and array accesses out of bounds
* non-finite math: typically due to divisions by zero, square roots of
  negative numbers or logarithms of negative numbers
* unhandled exception: always fatal when running with multiple MPI ranks

Many algorithms require parameters to be provided within valid ranges.
Range checks are implemented to catch invalid input values and generate
meaningful error messages, however these checks cannot always catch errors
arising from an invalid combination of two or more features. If you encounter
issues with a script, you can activate extra runtime checks by enabling C++
assertions. This is achieved by updating the CMake project and rebuilding
|es| with:

.. code-block:: bash

    cmake . -D CMAKE_BUILD_TYPE=RelWithAssert
    make -j$(nproc)
    ./pypresso script.py

The resulting build will run slightly slower, but will produce an error
message for common issues, such as divisions by zero, array access out
of bounds, or square roots of negative numbers.

If this still doesn't help, activate debug symbols to help with instrumentation:

.. code-block:: bash

    cmake . -D CMAKE_BUILD_TYPE=Debug
    make -j$(nproc)
    ./pypresso script.py 2>&1 | c++filt

The resulting build will be quite slow but segmentation faults will generate
a complete backtrace, which can be parsed by ``c++filt`` to demangle symbol
names. If this is not sufficient to track down the source of the error,
a debugging tool like GDB can be attached to |es| to catch the segmentation
fault signal and generate a backtrace. See :ref:`using GDB<GDB>` for more details.

If you are dealing with a segmentation fault or undefined behavior, and GDB
doesn't help or is too cumbersome to use (e.g. in MPI-parallel simulations),
you can as a last resort activate sanitizers:

.. code-block:: bash

    cmake . -D ESPRESSO_BUILD_WITH_ASAN=ON \
            -D ESPRESSO_BUILD_WITH_UBSAN=ON \
            -D CMAKE_BUILD_TYPE=RelWithAssert
    make -j$(nproc)
    ./pypresso script.py

The resulting build will be around 5 times slower that a debug build,
but it will generate valuable reports when detecting fatal exceptions.

It is possible to attach an external debugger to ``pypresso``, albeit with
a custom syntax. The ``pypresso`` executable file is actually not a program
but a script which sets the Python path appropriately and starts the Python
interpreter with user-defined arguments. Thus it is not possible to directly
run ``pypresso`` in a debugger; instead one has to use pre-defined command
line options:

.. code-block:: bash

     ./pypresso --tool script.py

where ``--tool`` can be any tool from the :ref:`table below <Debugging es with tools>`.
Only one tool can be used at a time. Some tools benefit from specific build
options, as outlined in the sections that follow. Most tools accept arguments
``<args>`` via the following variant:

.. code-block:: bash

     ./pypresso --tool="<args>" script.py

The sequence or arguments is passed as a string, which will be split at
whitespace characters by the shell interpreter. When the arguments need
whitespaces or quotation marks, those need to be properly escaped. When
no arguments are passed, sensible default values will be used instead.

.. _Debugging es with tools:

.. table:: Tools for the Python wrapper to |es|.

    +------------------------+-------------------------------------------------------------+
    | Tool                   | Effect                                                      |
    +========================+=============================================================+
    | ``--gdb``              | ``gdb --args python script.py``                             |
    +------------------------+-------------------------------------------------------------+
    | ``--lldb``             | ``lldb -- python script.py``                                |
    +------------------------+-------------------------------------------------------------+
    | ``--valgrind``         | ``valgrind --leak-check=full python script.py``             |
    +------------------------+-------------------------------------------------------------+
    | ``--cuda-gdb``         | ``cuda-gdb --args python script.py``                        |
    +------------------------+-------------------------------------------------------------+
    | ``--cuda-memcheck``    | ``cuda-memcheck python script.py``                          |
    +------------------------+-------------------------------------------------------------+
    | ``--cuda-sanitizer``   | ``compute-sanitizer --leak-check full python script.py``    |
    +------------------------+-------------------------------------------------------------+
    | ``--kernprof``         | ``kernprof --line-by-line --view script.py``                |
    +------------------------+-------------------------------------------------------------+

.. _Profiling:

Profiling
~~~~~~~~~

|es| is designed to leverage highly parallel computing environments and GPU
accelerators. To facilitate the investigation of communication bottlenecks
and inefficient algorithms, several profilers are natively supported,
with annotation markers placed in performance-critical parts of the C++ core.

.. _GDB:

GDB
~~~

.. note::

    Requires a debug build, enabled with the CMake option
    ``-D CMAKE_BUILD_TYPE=Debug``, as well as an external dependency:

    .. code-block:: bash

        sudo apt install gdb

The GNU Debugger (GDB) :cite:`stallman11a` is used to observe and control
the execution of C++ applications. GDB can catch signals, suspend the
program execution at user-defined break points, expose the content of
C++ variables and run C++ functions that have no side effects.

Here is a typical GDB session. Runs the failing simulation
with the pypresso ``--gdb`` flag to attach the process to GDB.
To catch a runtime error, use e.g. ``catch throw std::runtime_error``.
To catch a specific function, use ``break`` followed by the function name
(answer yes to the prompt about pending the breakpoint), or alternatively
provide the absolute filepath and line number separated by a colon symbol.
For a segmentation fault, no action is needed since it is automatically
caught via the SIGSEV signal. Run the simulation with ``run`` and wait
for GDB to suspend the program execution. At this point, use ``bt`` to
show the complete backtrace, then use ``frame <n>`` with ``<n>`` the number
of the innermost frame that is located inside the |es| source directory,
and finally use ``tui e`` to show the offending line in the source code
(``tui d`` to hide the source code). Use ``up`` and ``down`` to move in
the backtrace. The value of local variables can be inspected by GDB.
For a self-contained example, see the :ref:`GDB example<GDB-example>`.

It is possible to debug an MPI-parallel simulation script with GDB.
Keep in mind that contrary to a textbook example MPI application, where
all ranks execute the ``main`` function, in |es| the worker nodes are idle
until the head node on MPI rank 0 delegates work to them. This means that
on MPI rank > 1, break points will only have an effect in code that can be
reached from a callback function whose pointer has been registered in the
:ref:`MPI callbacks framework <The MPI callbacks framework>`.

The following command runs a script with 2 MPI ranks and binds a terminal
to each rank:

.. code-block:: bash

    mpiexec -np 2 xterm -fa 'Monospace' -fs 12 -e ./pypresso --gdb simulation.py

It can also be done via ssh with X-window forwarding:

.. code-block:: bash

    ssh -X username@hostname
    mpiexec -n 2 -x DISPLAY="${DISPLAY}" xterm -fa 'Monospace' -fs 12 \
        -e ./pypresso --gdb simulation.py

The same syntax is used for C++ unit tests:

.. code-block:: bash

    mpiexec -np 2 xterm -fa 'Monospace' -fs 12 \
        -e gdb src/core/unit_tests/EspressoSystemStandAlone_test

.. _GDB-example:

**GDB example**

To recreate a typical debugging session, let's purposefully introduce a null
pointer dereference in the ``int integrate()`` function, like so:

.. code-block:: c++

    int integrate(int n_steps, int reuse_forces) {
      int test = *std::shared_ptr<int>();

Running any simulation should produce the following trace:

.. code-block:: none

    $ ./pypresso ../samples/lj_liquid.py 2>&1 | c++filt
    *** Process received signal ***
    Signal: Segmentation fault (11)
    Signal code: Address not mapped (1)
    Failing at address: (nil)
    [ 0] /lib/x86_64-linux-gnu/libc.so.6(+0x42520)
    [ 1] /home/user/espresso/build/src/core/espresso_core.so(integrate(int, int)+0x49)
    [ 2] /home/user/espresso/build/src/core/espresso_core.so(integrate_with_signal_handler(int, int, bool)+0xaf)

Running in GDB should automatically catch the SIGSEV signal and allow us to
inspect the code and the state of all local variables:

.. code-block:: none

    $ ./pypresso --gdb ../samples/lj_liquid.py
    (gdb) run
    Thread 1 "python3.10" received signal SIGSEGV, Segmentation fault.
    in integrate (n_steps=20, reuse_forces=-1)
    at /home/user/espresso/src/core/integrate.cpp:260
    260   int test = *std::shared_ptr<int>();
    (gdb) bt
    #0  in integrate (n_steps=20, reuse_forces=-1)
        at /home/user/espresso/src/core/integrate.cpp:260
    #1  in integrate_with_signal_handler (n_steps=20, reuse_forces=-1,
          update_accumulators=false)
        at /home/user/espresso/src/core/integrate.cpp:484
    #2  in ScriptInterface::Integrators::SteepestDescent::integrate (
          this=..., params=std::unordered_map with 1 element = {...})
        at /home/user/espresso/src/script_interface/integrators/SteepestDescent.cpp:44
    (gdb) frame 0
    #0  in integrate (n_steps=20, reuse_forces=-1)
        at /home/user/espresso/src/core/integrate.cpp:260
    260   int test = *std::shared_ptr<int>();
    (gdb) tui e
    ┌─/home/user/espresso/src/core/integrate.cpp───────────────────────────────────┐
    │      257  }                                                                  │
    │      258                                                                     │
    │      259  int integrate(int n_steps, int reuse_forces) {                     │
    │  >   260    int test = *std::shared_ptr<int>();                              │
    │      261                                                                     │
    │      262    // Prepare particle structure and run sanity checks              │
    │      263    on_integration_start(time_step);                                 │
    └──────────────────────────────────────────────────────────────────────────────┘
    (gdb) print n_steps
    $1 = 20
    (gdb) ptype time_step
    type = double

.. _CUDA_GDB:

CUDA-GDB
~~~~~~~~

.. note::

    Requires a CUDA debug build, enabled with the CMake options
    ``-D ESPRESSO_BUILD_WITH_CUDA=ON -D CMAKE_BUILD_TYPE=Debug``.

The CUDA-GDB debugger :cite:`misc-cuda-gdb` is used to observe and control
the execution of CUDA applications. CUDA-GDB can catch signals, suspend the
program execution at user-defined break points and expose values in CUDA
variables. When a signal is caught inside a CUDA kernel, the stack trace
only shows device function calls. When stepping into a CUDA kernel launch,
the stack trace shows both host and device function calls.

.. _ASAN:

ASAN
~~~~

.. note::

    Requires specific compiler and linker flags, enabled with the CMake option
    ``-D ESPRESSO_BUILD_WITH_ASAN=ON -D CMAKE_BUILD_TYPE=RelWithAssert``.

The AddressSanitizer (ASAN) :cite:`serebryany12a` is a memory error detection
tool. It detects memory leaks and bugs caused by dangling references.

For more details, please consult the tool online documentation [5]_.

.. _UBSAN:

UBSAN
~~~~~

.. note::

    Requires specific compiler and linker flags, enabled with the CMake option
    ``-D ESPRESSO_BUILD_WITH_UBSAN=ON -D CMAKE_BUILD_TYPE=RelWithAssert``.

The UndefinedBehaviorSanitizer (UBSAN) :cite:`misc-ubsan` is a detection tool
for undefined behavior. It detects bugs caused by dangling references,
array accesses out of bounds, signed integer overflows, etc.

For more details, please consult the tool online documentation [6]_.

.. _Caliper:

Caliper
~~~~~~~

.. note::

    Requires external features ``CALIPER``, enabled with the CMake option
    ``-D ESPRESSO_BUILD_WITH_CALIPER=ON``.

Caliper [1]_ :cite:`boehme16a` is a low-overhead annotation library for C++.
By default, |es| comes with several markers in performance-critical parts
of the main integration loop.

In the example below, a P3M simulation is profiled to reveal that the
short-range loop (N-squared summation for Lennard-Jones and Coulomb)
and long-range forces (FFT summation) contribute equally to the runtime:

.. code-block:: none

    $ CALI_CONFIG_PROFILE=runtime-report ./pypresso ../samples/p3m.py --cpu
    Path                         Inclusive time Exclusive time    Time %
    integrate                             14.18           0.01      0.08
      Integration loop                    13.84           0.43      2.88
        force_calc                        13.41           0.20      1.35
          copy_forces_from_GPU             0.01           0.01      0.07
          short_range_loop                 6.55           6.55     44.02
          calc_long_range_forces           6.40           6.40     43.00
          init_forces                      0.24           0.24      1.58
          copy_particles_to_GPU            0.01           0.01      0.07

For the GPU implementation of the P3M algorithm, the long-range force
calculation is cheaper, however the transfer of particle data to and from
the GPU incur additional costs that are not negligible:

.. code-block:: none

    $ CALI_CONFIG_PROFILE=runtime-report ./pypresso ../samples/p3m.py --gpu
    Path                         Inclusive time Exclusive time    Time %
    integrate                             14.30           0.03      0.14
      Integration loop                    13.87           1.76      7.90
        force_calc                        12.12           0.82      3.68
          copy_forces_from_GPU             2.09           2.09      9.42
          short_range_loop                 3.20           3.20     14.38
          calc_long_range_forces           3.75           3.75     16.87
          init_forces                      1.25           1.25      5.61
          copy_particles_to_GPU            1.01           1.01      4.56

For a more fine-grained report on GPU kernels:

.. code-block:: none

    $ CALI_CONFIG=cuda-activity-report ./pypresso ../samples/p3m.py --gpu

To introduce custom markers at the C++ level, add ``CALI`` macros inside
performance-critical functions to register them:

.. code-block:: c++

    void force_calculation(CellStructure &cell_structure, double time_step) {
    #ifdef CALIPER
      CALI_CXX_MARK_FUNCTION;
    #endif
      /* ... */
    }

To introduce custom markers at the Python level,
use a :class:`~espressomd.profiler.Caliper` object to fence code blocks:

.. code-block:: python

    import espressomd.profiler
    cali = espressomd.profiler.Caliper()
    cali.begin_section(label="calc_energies")
    energies = system.analysis.energy()
    cali.end_section(label="calc_energies")

.. _Valgrind:

Valgrind
~~~~~~~~

.. note::

    Requires external features ``VALGRIND`` and debug symbols,
    enabled with the CMake options
    ``-D ESPRESSO_BUILD_WITH_VALGRIND=ON -D CMAKE_BUILD_TYPE=RelWithDebInfo``,
    as well as external dependencies:

    .. code-block:: bash

        sudo apt install valgrind kcachegrind graphviz
        python3 -m pip install --user gprof2dot

The Valgrind [2]_ :cite:`nethercote07a,nethercote03a` framework brings several
tools to examine a program runtime performance.

.. _Callgrind:

Callgrind
"""""""""

The Callgrind [3]_ :cite:`weidendorfer04a` tool generates a graph of function
calls. This type of instrumentation has a lot of overhead, therefore the time
spent in functions might not always be reliable, and the program execution
is slowed down significantly. To remediate the latter, it is common to
restrict instrumentation to a specific part of the code using markers.
By default, |es| comes with markers in the integration loop,
which is the most performance-critical part of the core.

In the following example, the P3M algorithm is profiled to generate a call
graph that can be converted to a static graph using ``gprof2dot`` and ``dot``:

.. code-block:: bash

    ./pypresso --valgrind="--tool=callgrind --instr-atstart=no" ../samples/p3m.py --cpu
    callgrind_out=$(ls -t -1 callgrind.out.*[[:digit:]] | head -1)
    python3 -m gprof2dot --format=callgrind --output=${callgrind_out}.dot ${callgrind_out}
    dot -Tpdf ${callgrind_out}.dot -o ${callgrind_out}.pdf

The Valgrind output file generally follows the pattern ``callgrind.out.pid``,
where ``pid`` is the actualy process id. The ``${callgrind_out}`` variable
is populated with the return value of a subshell commands that finds the most
recent output file that matches that pattern.

It is also possible to open the output file in KCachegrind [4]_ to browse
the call graph interactively and visualize the time spent in each function:

.. code-block:: bash

    kcachegrind ${callgrind_out}

.. _Compute Sanitizer:

Compute Sanitizer
~~~~~~~~~~~~~~~~~

.. note::

    Requires a CUDA build, enabled with the CMake options
    ``-D ESPRESSO_BUILD_WITH_CUDA=ON``.

The Compute Sanitizer [9]_ :cite:`misc-compute-sanitizer` framework is similar
to :ref:`Valgrind`, but for NVIDIA GPUs. The exact command line options
differ with the CUDA version. If the command line examples below don't work,
please refer to the NVIDIA user guide version that corresponds to the locally
installed CUDA toolkit.

To detect memory leaks:

.. code-block:: bash

    ./pypresso --cuda-sanitizer="--tool memcheck --leak-check full" script.py

Add option ``--error-exitcode 1`` to return an error code when issues are detected.

To detect access to uninitialized data:

.. code-block:: bash

    ./pypresso --cuda-sanitizer="--tool initcheck" script.py

Checking for uninitialized data is quite expensive
for the GPU and can slow down other running GPU processes.

.. _perf:

perf
~~~~

.. note::

    Requires debug symbols, enabled with the CMake option
    ``-D CMAKE_BUILD_TYPE=DebugOptimized``,
    as well as external dependencies. On Ubuntu:

    .. code-block:: bash

        sudo apt install linux-tools-generic

    On Debian:

    .. code-block:: bash

        sudo apt install linux-perf

    On Fedora:

    .. code-block:: bash

        sudo apt install perf

The perf [7]_ :cite:`misc-perf` tool generates a graph of function calls
with time measurements.
It requires privileges that can only be set as root.

In the following example, the P3M algorithm is profiled to generate a call
graph in a file called ``perf.data``, which is then read to generate a report:

.. code-block:: bash

    original_value=$(sysctl -n kernel.perf_event_paranoid)
    sudo sysctl -w kernel.perf_event_paranoid=3
    perf record --call-graph dwarf ./pypresso ../samples/p3m.py --cpu
    sudo sysctl -w kernel.perf_event_paranoid=${original_value}
    perf report --call-graph

When inside the report, press ``/`` to search for a function name,
e.g. ``integrate``, then highlight the symbol and press ``+`` to expand
its call graph. Press ``q`` to exit the program, or close open tabs.

A large amount of data will be written to disk during the recording step,
typically several hundred megabytes. If the hard drive write latency
is too high, the following warning will be emitted:

.. code-block:: none

    Warning:
    Processed 17655 events and lost 7 chunks!
    Check IO/CPU overload!

Using a tmpfs drive, perf can write the file directly to RAM
(mounted as a filesystem), which has better latency.
To get a list of mounted tmpfs drives and their capacity:

.. code-block:: none

    $ mount | grep "tmpfs"
    tmpfs on /dev/shm type tmpfs (rw,nosuid,nodev)
    $ df -h /dev/shm/
    Filesystem      Size  Used Avail Use% Mounted on
    tmpfs            32G  320K   32G   1% /dev/shm

To use a tmpfs drive as storage:

.. code-block:: bash

    perf record --call-graph dwarf -o /dev/shm/perf.data ../samples/p3m.py --cpu
    perf report --call-graph -i /dev/shm/perf.data
    rm /dev/shm/perf.data

.. _kernprof:

kernprof
~~~~~~~~

.. note::

    Requires an external dependency:

    .. code-block:: bash

        python3 -m pip install --user line_profiler

kernprof [8]_ :cite:`misc-kernprof` is a low-overhead Python profiler.
It supports two instrumentation modes: ``line_profile`` and ``cProfile``.
The ``--builtin`` option injects a ``LineProfiler`` object and a ``profile``
function in the global namespace of the instrumented script.
The latter can be used as a decorator (``@profile``),
as a context manager (``with profile:``), or
as begin/end markers (``profile.enable()``, ``profile.disable()``)
to select the regions of code to instrument,
although the ``line_profile`` mode only supports the decorator behavior.
The ``line_profile`` mode cannot instrument code from imported modules,
whereas the ``cProfile`` mode can.

To make the instrumented script executable with and without kernprof
when using decorators, add the following code at the top of the script:

.. code-block:: python

    if "line_profiler" not in dir():
        def profile(func):
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper

To run kernprof in ``line_profile`` mode:

.. code-block:: bash

    ./pypresso --kernprof="--line-by-line --view" ../samples/p3m.py --cpu

To later view the results again:

.. code-block:: bash

    python3 -m line_profiler p3m.py.lprof

To run kernprof in ``cProfile`` mode:

.. code-block:: bash

    ./pypresso --kernprof="" ../samples/p3m.py --cpu

To interactively read the data:

.. code-block:: none

    python3 -m pstats p3m.py.prof
    p3m.py.prof% sort time
    p3m.py.prof% reverse
    p3m.py.prof% stats
      ncalls  tottime  percall  cumtime  percall filename:lineno(function)
           2  1.090    0.545    1.090    0.545   /opt/espressomd/integrate.py:156(run)
           1  1.817    1.817    1.817    1.817   /opt/espressomd/electrostatics.py:71(_activate)
          10  2.619    0.262    2.619    0.262   /opt/espressomd/integrate.py:101(run)
    p3m.py.prof% quit

____

.. [1]
   https://software.llnl.gov/Caliper/

.. [2]
   https://valgrind.org/docs/manual/

.. [3]
   https://valgrind.org/docs/manual/cl-manual.html

.. [4]
   https://kcachegrind.github.io/html/Home.html

.. [5]
   https://github.com/google/sanitizers/wiki/AddressSanitizer

.. [6]
   https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html

.. [7]
   https://perf.wiki.kernel.org/index.php/Main_Page

.. [8]
   https://github.com/pyutils/line_profiler

.. [9]
   https://docs.nvidia.com/compute-sanitizer/ComputeSanitizer/index.html
