
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


.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
