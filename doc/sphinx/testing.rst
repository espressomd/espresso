
.. _Unit testing:

Unit testing
============

Tests are configured in CMake and run by CTest.


.. _Cpp unit tests:

C++ unit tests
--------------

Testing individual C++ functions and classes.


Running tests
^^^^^^^^^^^^^

To build and run the tests:

.. code-block:: bash

   make -j$(nproc) check_unit_tests


Writing tests
^^^^^^^^^^^^^

Framework: `Boost.Test <https://www.boost.org/doc/libs/release/libs/test/>`_.

Store new unit tests in :file:`/src/core/unit_tests/`. Tests usually follow
this structure:

.. code-block:: c++

   /// @file /src/core/unit_tests/Particle_test.cpp
   #include <boost/test/unit_test.hpp>
   #include "core/utils/serialization/Particle.hpp"

   BOOST_AUTO_TEST_CASE(comparison) {
     {
       Particle p, q;
       p.identity() = 1;
       q.identity() = 2;
       BOOST_CHECK(p != q);
     }
     {
       Particle p, q;
       p.identity() = 2;
       q.identity() = 2;
       BOOST_CHECK(p == q);
     }
   }


Mocked classes should be enclosed in a namespace and stored in
:file:`/src/core/unit_tests/mock`. Here is an example of a mocked Cell:

.. code-block:: c++

   /// @file /src/core/unit_tests/mock/Cell.hpp
   namespace Testing {
     template <typename Particle> class Cell {
       public:
       Cell() : n(0) {}
       std::size_t n;
       std::vector<Particle> part;
     };
   } // namespace Testing

which is then used as:

.. code-block:: c++

   #include "mock/Cell.hpp"
   using Cell = Testing::Cell<Particle>;


CMake
^^^^^

Configure tests in :file:`/src/core/unit_tests/CMakeLists.txt` using the syntax

.. code-block:: cmake

   unit_test(NAME <mytest> SRC <mytest.cpp> [DEPENDS <target1>[, ...]] [NUM_PROC <N>])

where ``NUM_PROC`` instructs CTest to run the binary through OpenMPI with ``N``
threads. When ``NUM_PROC`` is not provided, the binary is executed normally.
The use of more than 1 source file should be avoided in favor of the
``DEPENDS`` option.


.. _Python integration tests:

Python integration tests
------------------------

Testing Python bindings and numerical results of core features.


Running tests
^^^^^^^^^^^^^

To build and run the tests:

.. code-block:: bash

   make -j$(nproc) check_python_serial
   make -j$(nproc) check_python_parallel
   make -j$(nproc) check_python_parallel_odd


Writing tests
^^^^^^^^^^^^^

Framework: `unittest <https://docs.python.org/3/library/unittest.html>`_.

Store new unit tests in :file:`/testsuite/python/`. Tests usually follow
this structure:

.. code-block:: python

   # /testsuite/python/constraint_shape_based.py
   import unittest as ut
   import numpy as np
   import espressomd

   @ut.skipIf(not espressomd.has_features(["LENNARD_JONES_GENERIC"]),
              "Features not available, skipping test!")
   class ShapeBasedConstraintTest(ut.TestCase):

       box_l = 30.
       system = espressomd.System(box_l=3 * [box_l])

       def tearDown(self):
           self.system.part.clear()
           self.system.constraints.clear()

       def test_hollowcone(self):
           system = self.system
           system.time_step = 0.01
           system.cell_system.skin = 0.4
           <...>

   if __name__ == "__main__":
       ut.main()


Python decorators ``@ut.skipIf(<condition>)`` are used to skip classes and
methods when optional features are not compiled in.


CMake
^^^^^

Configure tests in :file:`/testsuite/python/CMakeLists.txt` using the syntax

.. code-block:: cmake

   python_test(FILE <mytest.py>
               MAX_NUM_PROC <N>
               [DEPENDENCIES <../dependency1.py>[, ...]]
               [CONFIGURATIONS <configuration_list>])

where ``MAX_NUM_PROC`` instructs CTest to run the script through OpenMPI
with at most ``N`` threads. The actual number of threads used depends on
``CONFIGURATIONS``: ``"serial"`` (one core), ``"parallel"`` (even number
of cores) and ``"parallel_odd"`` (odd number of cores). The actual value
is determined by CMake. Files listed in ``DEPENDENCIES`` are passed to
``configure_file()``.


.. _IPython notebooks and Python samples:

IPython notebooks and Python samples
------------------------------------

Introduction
^^^^^^^^^^^^

:file:`/samples/` contains Python scripts showcasing typical uses of |es|.
:file:`/doc/tutorials/` contains IPython notebooks and bonus Python scripts
used in teaching sessions. They are both part of a different CI schedule.


Running tests
^^^^^^^^^^^^^

To run the tests in parallel:

.. code-block:: bash

   make -j$(nproc) check_tutorials ARGS=-j$(nproc)
   make -j$(nproc) check_samples ARGS=-j$(nproc)


Writing tests
^^^^^^^^^^^^^

Tutorials and samples are designed for interactive use, and cannot be
imported like conventional Python modules. For example, IPython notebooks have
to be first converted to Python scripts, then imported with the ``importlib``
module because of their non-standard filenames. Some scripts need to be
imported together in the same Python session for them to work, while others
need to access resources (i.e. .dat files) found in the same directory as
the scripts. This last issue is solved by copying the complete
:file:`/samples` directory to
:file:`/build/testsuite/scripts/samples/local_samples` and the
:file:`/doc/tutorials` directory to
:file:`/build/testsuite/scripts/tutorials/local_tutorials`.

Since importing a Python script causes its execution, the simulation will
run upon import, at which point all global variables become accessible to
the unittest classes, including the :class:`~espressomd.system.System` object.
However, some scripts can be slow to import and require specific command line
arguments. To solve both issues, the scripts are edited first and saved to
a new file with suffix ``_processed.py``, which is the one being actually
imported by the testing script. During editing, global variables controlling
the running time (number of integration steps, number of particles, target
accuracy, etc.) are substituted with new values defined in the testing script,
and ``sys.argv`` is modified to contain the command line arguments. Several GUI
classes are also disabled during this step.

Here is a test template to load ``sample.py``, while altering the values of two
global variables (``warm_steps``, ``n_iterations=20``), setting up command line
arguments (``--cpu 0.001``) and replacing an infinite loop with a finite loop:

.. code-block:: python

   import unittest as ut
   import importlib_wrapper

   def disable_visualizer_GUI(code):
       breakpoint = "while True:"
       assert breakpoint in code
       code = code.replace(breakpoint, "for _ in range(5):", 1)

   sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
       "@SAMPLES_DIR@/sample.py", cmd_arguments=["--cpu", "0.001"],
       warm_steps=100, n_iterations=20, substitutions=disable_visualizer_GUI)

   @skipIfMissingFeatures
   class Sample(ut.TestCase):
       system = sample.system

       def test_something(self):
           self.assertLess(abs(sample.pressure - 1.0), 1e-3)
           <...>

   if __name__ == "__main__":
       ut.main()


Contrary to the :ref:`Python integration tests` where fixtures
``@ut.skipIf(<condition>)`` are used to disable tests for missing features,
tutorials and samples already have ``espressomd.assert_features()`` statements
from which a fixture ``@skipIfMissingFeatures`` is automatically created.
When importing multiple scripts in the same test (for example, some tutorials
have bonus scripts that don't necessarily require ``assert_features()``),
on the first import failure all subsequent imports will be skipped.

Please note that numerical results of interest (``sample.pressure`` in the
previous example) need to be stored in global variables to be accessible.
It is also important to format optional code cells in IPython/Jupyter
notebooks as **markdown cells** to prevent them from running during import.
Simply enclose them in triple backticks to preserve syntax highlighting:

.. code-block:: md

   Alternatively, you could start the visualizer using:

   ```python
   from espressomd import visualization
   visualizer = visualization.openGLLive(system)
   visualizer.run()
   ```


Given the stochastic nature of the scripts (NumPy and |es| RNG seeds are
assigned randomly), it is necessary to run new tests multiple times to ensure
reproducibility of the numerical results being tested. Here is one method:

.. code-block:: bash

   for i in {1..50}
   do
     ../../../pypresso test_drude_bmimpf6_with_cpu.py
     if [ $? != "0" ]
     then
       echo $i
       break
     fi
   done


Store new tests in :file:`/testsuite/scripts/samples/` or
:file:`/testsuite/scripts/tutorials/`.


CMake
^^^^^

Configure tests in :file:`/testsuite/scripts/samples/CMakeLists.txt` or
:file:`/testsuite/scripts/tutorials/CMakeLists.txt` using the syntax

.. code-block:: cmake

   sample_test(FILE <mytest.py>
               [LABELS "<label1>[; ...]"]
               [SUFFIX "<suffix>"]
               [DEPENDENCIES <../dependency1.py>[, ...]])

for samples, or ``tutorial_test()`` for tutorials and their bonus scripts.

Files listed in ``DEPENDENCIES`` are passed to ``configure_file()``.

The ``GPU`` flag instructs CTest to skip the test if no CUDA-capable GPU is
available, or to run GPU tests one after the other if CTest runs in parallel
(by default GPU rank 0 is used, so we can only run one GPU job at a time).
This example shows how to set up an LB test to use the CPU implementation once
and the GPU implementation once (if a GPU is available):

.. code-block:: cmake

   sample_test(FILE test_lbf.py SUFFIX cpu)
   sample_test(FILE test_lbf.py SUFFIX gpu LABELS "gpu")

where ``SUFFIX`` is used to create a file :file:`test_lbf_with_cpu.py` and a
file :file:`test_lbf_with_cpu`. In general, ``SUFFIX`` is used to generate
multiple test scripts from a template. This example shows how to generate tests
for various :class:`~espressomd.shapes` using a loop:

.. code-block:: cmake

   foreach(shape wall;sphere;ellipsoid;cylinder;hollowcone)
     sample_test(FILE test_visualization_constraints.py SUFFIX ${shape})
   endforeach(shape)


Sequentiality can be enforced with fixtures:

.. code-block:: cmake

   sample_test(FILE test_save_checkpoint.py)
   sample_test(FILE test_load_checkpoint.py)
   set_tests_properties(sample_save_checkpoint PROPERTIES FIXTURES_SETUP    saved_checkpoint)
   set_tests_properties(sample_load_checkpoint PROPERTIES FIXTURES_REQUIRED saved_checkpoint)


.. _Installation tests:

Installation tests
------------------

Introduction
^^^^^^^^^^^^

Test the installation of |es| and its Python bindings.


Running tests
^^^^^^^^^^^^^

To run the tests:

.. code-block:: bash

   make check_cmake_install


Writing tests
^^^^^^^^^^^^^

Framework: custom Bash script :file:`/testsuite/cmake/BashUnitTests.sh`.

Here is a toy example to check if a file exists and if a Python module can be
imported:

.. code-block:: bash

   #!/usr/bin/env bash

   # load bash unit testing library
   source BashUnitTests.sh

   function test_install() {
     assert_file_exists "@CMAKE_BINARY_DIR@/ipypresso"
   }

   function test_import() {
     local import_dir="@DESTDIR@/@CMAKE_INSTALL_PREFIX@/@Python_SITEARCH@"
     local instruction="import sys;sys.path.insert(0, '${import_dir}');import espressomd"
     assert_return_code "@CMAKE_BINARY_DIR@/pypresso" -c "${instruction}"
   }

   # run tests
   run_test_suite

Store new tests in :file:`/testsuite/cmake/` and add Execute permission with
``chmod +x test_<name>.sh`` to avoid the following CTest error message:

.. code-block:: none

   The following tests FAILED:
         1 - test_python_bindings (BAD_COMMAND)


CMake
^^^^^

Configure tests in :file:`/testsuite/cmake/CMakeLists.txt` using the syntax

.. code-block:: cmake

   cmake_test(FILE <mytest.sh> [DEPENDENCIES <../dependency1.sh>[, ...]])

Files listed in ``DEPENDENCIES`` are passed to ``configure_file()``.

