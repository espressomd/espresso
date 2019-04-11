
.. _Contact the Developers:

Contact the Developers
======================

To contact the |es| developers, please write an email to the developers mailing list:
espressomd-devel@nongnu.org
to subscribe to the developers' mailing list go to
http://lists.nongnu.org/mailman/listinfo/espressomd-devel


.. _Before you start a development project:

Before you start a development project
--------------------------------------
Before you start a development project for |es|, please always write to the developers mailing list and describe the project.
This is to avoid that several people work on the same thing at the same time. Also, implementation details can be discussed in advance. In many cases, existing developers can point to re-usable code and simpler solutions.


.. _Development Environment:

Development Environment
=======================


.. _Required Development Tools:

Required Development Tools
--------------------------

-  First of all, please install the dependencies for compiling |es|. See the section on "Getting, compiling and running" in the user guide.

-  To be able to access the development version of |es|, you will need
   the distributed versioning control system git_.

-  To build the sphinx documentation, you will need the Python packages listed in :file:`requirements.txt` in the top-level source directory. To install them, issue:

   .. code-block:: bash

      pip install --upgrade --user -r requirements.txt

   Note, that some distributions now use ``pip`` for Python3 and ``pip2`` for Python 2.

-  To build the tutorials, you will need LaTeX.

-  To compile the Doxygen code documentation, you will need to have the
   tool doxygen_.

All of these tools should be easy to install on most Unix operating
systems.

.. _Getting the Development Code:

Getting the Development Code
----------------------------
We use Github for storing the source code and its history, and for managing the development process.
The repository is located at http://github.com/espressomd/espresso.
To get the current development code, run:

.. code-block:: bash

  git clone git://github.com/espressomd/espresso

This will create a directory named "espresso" which contains the code.
The build process does not differ from the one for release versions described in the users' guide.


Build System
------------

The build system of |es| is based on CMake.

The central source files of the build system are the following:

-  :file:`CMakeLists.txt`

-  Contents of the :file:`cmake` directory

-  The :file:`CMakeLists.txt` files in the :file:`src/`, :file:`doc/`, and :file:`testsuite/` directories and their sub-directories

The most common reasons for editing these files are:

-  Adding new source files
-  Adding new external dependencies

Adding New Source Files
~~~~~~~~~~~~~~~~~~~~~~~

To add new files to |es| (like C++ source files or header files) you
need to look at the CMakeLists.txt in the directory where the file is located.

* Please note that .hpp-header files usually do not have to be added to CMakeList.txt

* All files are explicitly included (e.g., testsuite/python/CMakeLists.txt):: 

      set(py_tests  bondedInteractions.py
                   cellsystem.py
                   constraint_shape_based.py
                   coulomb_cloud_wall.py)


Testsuite
---------

-  New or significantly changed features will only be accepted, if they have a test case.
   This is to make sure, the feature is not broken by future changes to |es|, and so other users can get an impression of what behavior is guaranteed to work.
-  There are two kinds of tests:

  -  C++-unit tests, testing individual C++ functions and classes. They make use of the boost unit test framework and reside in :file:`src/core/unit_tests`
  -  Python integration tests, testing the Python interface and (physical) results of features. They reside in :file:`python`

- To execute the tests, run:

  .. code-block:: bash

     make check

  in the top build directory.


.. _Documentation:

Documentation
=============

The documentation of |es| consists of four parts:

  -  The users' guide and developers' guide are located in :file:`doc/sphinx`, and make use of the Sphinx Python package
  -  In-code documentation for the Python interface is located in the various files in src/python/espressomd and also makes use of the Sphinx Python package. We make use of the napoleon extension and use the NumPy documentation style.
  -  In-code documentation of the C++ core is located in the .cpp and .hpp files in :file:`/src/core` and its sub-directories and makes use of Doxygen.

Doxygen Code Documentation
--------------------------

The documentation of each function should contain a short description,
if necessary a more detailed description and a description for the
return value and parameters.

Look at the documentation of existing files and functions to get a
feeling how it should be!

Doxygen is able to understand simple LaTeX and HTML commands as well as
some special command in order to give the documentation a nice structure
and to make it more readable. In the following list you find a short
description of the most common commands we need:

-  | ``\anchor`` *name* *description*
   | Create an anchor to which you can refer using the ``\ref`` command.

-  | ``\ref`` *name* ``["``\ *text*\ ``"]``
   | Insert a link to another object in the documentation (*e.g.*\ an
     anchor).

-  | ``<a href="http://www.your_url.html">title</a>``
   | Link to an external HTML source.

-  | ``\file`` *name* *description*
   | Special anchor for a file.

-  | ``\image html`` *image*
   | Include a picture. The picture file should reside in the subdir
     ``doc/doxygen/figs``. Do not use the HTML ``<img>``-tag to include
     pictures, as doxygen_ will not copy the pictures into the
     documentation.

-  | ``<ul> <li>List entry 1</li> <li>List entry 2</li></ul>``
   | Creates a list in the documentation.

-  | ``\param`` *name* *description*
   | Document the parameter of a function.

-  | ``\return`` *description*
   | Document the return value of a function.

.. _Programmers's Guide:


Programmer's Guide
==================

This chapter provides some hints on how to extend |es|. It is not
exhaustive, so for major changes the best documentation are the other
developers.


Source code structure
---------------------

The source tree has the following structure:

* src: The actual source code

  * core: The C++ source code of the simulation core
  * python/espressomd: Source of the espressomd Python module and its submodules
  * script_interface: C++ source code of the script_interface component, which links Python classes to functionality in the simulation core

* doc: Documentation

  * sphinx: The sphinx-based documentation, consisting of user and developer guide.
  * tutorials/python: Source and pdf files for the introductory tutorials
  * doxygen_: Build directory for the C++ in-code documentation

* testsuite/python: Python integration tests. Note that some C++ unit tests for individual core components are in src/core/unittests
* samples/python: Some sample scripts
* libs: External dependencies (at this point h5xx)
* maintainer: Files used by the maintainers

  * configs: Collection of myconfig.hpp files which activate different sets of features for testing.
  * docker: Definitions of the docker images for various distributions used for continuous integration testing
  * CI: Support files for the continuous integration testing run on the Travis-CI service.
  * jenkins: Outdated support files for the Jenkins continuous integration testing


Flow control and communications architecture
--------------------------------------------
Espresso uses two communication models, namely master-slave and synchronous.

* When Espresso does not run an integration, it works in the master-slave mode, i.e. the head node (MPI rank 0) in a parallel simulation
  runs the Python script, whereas all other nodes are idle until they receive a command from the head node. Such commands include particle creation,
  changing of particle properties and changing global simulation parameters.
  When a Python command such as::

    system.part.add(pos=(1, 2, 3))

  is issued, the head node determines, which node is responsible for the given position, and then sends the node the command to place the particle.

* When an integration is started in Python on the head node, a command to start the integration is sent to all nodes, in the master-slave framework described above.
  Then, Espresso switches into the synchronous mode, in which all nodes run the same code in the integration loop at the same time.
  The code of the main integration loop is in ``integrate.cpp:integrate_vv()``.
  When writing code which is run during the main integration loop, no commands making use of the master-slave mechanism can be called.
  When code during the integration loop executes MPI communication, it has to be ensured, that the MPI call is executed on all nodes
  involved in the communication. If this is not done, a deadlock will result.

Adding calls to the master-slave framework
------------------------------------------

Using an instance of MpiCallback
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Write the callback slave function, which will be executed on all nodes except the head node (0):

  .. code-block:: c++

    void my_callback(int p1, int p2) {
      // Do something. The two int-parameters can be used for anything
    }

* On all nodes, the callback has to be registered:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add(my_callback);
    }

  You can, e.g., call your registration from ``initialize.cpp:on_program_start()``
  Instead of a static function, from which a ``std::function<void(int,int)>`` can be constructed can
  be used. For example:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void register_my_callback() {
      Communication::mpiCallbacks().add([](int, int){ /* Do something */ });
    }

  can be used to add a lambda function as callback.
* Then, you can use your callback from the head node:

  .. code-block:: c++

    #include "MpiCallbacks.hpp"
    void call_my_callback() {
      Communication::mpiCallbacks.call(my_callback, param1, param2);
    }

  This only works outside the integration loop. After the callback has been called, synchronous mpi communication can be done.

Legacy callbacks
~~~~~~~~~~~~~~~~

Older code uses callbacks defined in the ``CALLBACK_LIST`` preprocessor macro in :file:`communications.cpp`. They are called via ``mpi_call()``.
See ``communications.cpp:mpi_place_particle()`` for an example.

.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
