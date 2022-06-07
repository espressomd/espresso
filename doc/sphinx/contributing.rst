.. _Contributing:

Contributing
============

Up to date information about the development of |es| can be found in the
`ESPResSo wiki <https://github.com/espressomd/espresso/wiki>`__:

- proceedings of the ESPResSo meetings
- list of planned releases
- developer's guide

The official website at https://espressomd.org provides additional information:

-  Latest stable release of |es| and older releases
-  Obtaining development version of |es|
-  Archives of both developers' and users' mailing lists
-  Registering to mailing lists

.. _Contributing your own code:

Contributing your own code
--------------------------

If you are planning to make an extension to or already have a piece of
your own code which could be useful to others, you are very welcome to
contribute it to the community. Before you start making any changes to
the code, you should fork the `espressomd/espresso
<https://github.com/espressomd/espresso>`__ repository and work in a new branch.

It is also generally a good idea to contact the developers on the mailing lists
before you start major coding projects. It might be that someone else is already
working on the problem or has a solution at hand.

You will find more detailed information on our development processes in
`Contributing to ESPResSo
<https://github.com/espressomd/espresso/blob/python/CONTRIBUTING.md>`__.
Please also refer to our developer's guide in the
`ESPResSo wiki <https://github.com/espressomd/espresso/wiki>`__.

.. _Required Development Tools:

Required Development Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^

-  First of all, please install the dependencies for compiling |es|.
   See the section :ref:`Installation`.

-  To be able to access the development version of |es|, you will need the
   distributed versioning control system git_ and a GitHub account to fork the
   `espressomd/espresso <https://github.com/espressomd/espresso>`__ repository

-  To build the user documentation, you will need Sphinx_.

-  To build the tutorials, you will need Jupyter_.

-  To build the core documentation, you will need Doxygen_.

All of these tools should be easy to install on most Unix operating systems.

You can find all Python dependencies of |es| in :file:`requirements.txt` in the
top-level source directory. Several optional packages for graphics, external
devices and continuous integration (CI) are not strictly essential and can be
safely removed if you're planning on installing dependencies via ``pip``:

.. code-block:: bash

   pip3 install --upgrade --user -r requirements.txt

Note that some distributions now use ``pip3`` for Python3 and ``pip2`` for
Python2.


.. _Building the User's guide:

Building the User's guide
-------------------------

If, while reading this documentation, you notice any mistakes or
undocumented features, you are very welcome to
contribute to the guide and have others benefit from your knowledge.

For this, you should clone the development version at `espressomd/espresso
<https://github.com/espressomd/espresso>`__. Next build the software as shown
in :ref:`Installation` and then the documentation with ``make sphinx``.


.. _git: https://git-scm.com/

.. _Doxygen: https://www.doxygen.nl/manual/index.html

.. _Sphinx: https://www.sphinx-doc.org/en/master/

.. _Jupyter: https://jupyter.org
