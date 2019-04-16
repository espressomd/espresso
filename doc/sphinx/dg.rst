
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


Testsuite
---------

-  New or significantly changed features will only be accepted if they have a test case.
   This is to make sure the feature is not broken by future changes to |es|, and so other users can get an impression of what behavior is guaranteed to work.
-  There are multiple kinds of tests:

  -  C++-unit tests, testing individual C++ functions and classes. They make use of the boost unit test framework and reside in :file:`/src/core/unit_tests`
  -  Python integration tests, testing the Python interface and (physical) results of features. They reside in :file:`/testsuite/python`
  -  CMake tests, testing the software can be successfully installed. They reside in :file:`/testsuite/cmake`
  -  Python scripts tests, testing the IPython notebooks and Python samples in :file:`/doc/tutorials` and :file:`/samples`. They reside in :file:`/testsuite/scripts/tutorials` resp. :file:`/testsuite/scripts/samples` and are executed on a different schedule

- To execute the tests, run:

  .. code-block:: bash

     make check

  in the top build directory.

- See :ref:`Unit testing` for how to develop new tests

.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
