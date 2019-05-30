.. _Contributing:

Contributing
============

Up to date information about the development of |es| can be found at the web
page http://espressomd.org As the important information can change in
time, we will not describe its contents in detail but rather request the
reader to go directly to the URL. Among other things, one can find
information about the following topics there:

-  Latest stable release of and older releases

-  Obtaining development version of |es|

-  Archives of both developers' and users' mailing lists

-  Registering to mailing lists

-  Submitting a bug report

.. _Community support:

Community support
-----------------

If you have any questions concerning |es| which you cannot resolve by
yourself, you may search for an answer on:

-  the `GitHub issue tracker <https://github.com/espressomd/espresso/issues>`_
   (remove the ``is:open`` option in the search bar)
-  the `espressomd-users mailing list archive
   <http://lists.nongnu.org/archive/html/espressomd-users/>`_

If you still didn't find a solution, you may consider either opening a new issue
on GitHub or sending an email on the mailing list. Instructions on how to
register to the mailing lists and post messages can be found in `Mailing Lists
<http://espressomd.org/wordpress/community-and-support/mailing-lists/>`_.

For several reasons it is recommended to send all questions to the issue
tracker or mailing list rather than to contact individual developers:

-  All registered users get your message and you have a higher
   probability that it is answered soon.
-  Your question and the answers are archived and the archives can be
   searched by others.
-  The answer may be useful also to other registered users.
-  There may not be a unique answer to your problem and it may be useful
   to get suggestions from different people.

Please remember that this is a community mailing list and a community issue
tracker. It is other users and developers who are answering your questions.
They do it in their free time and are not paid for doing it.

.. _Contributing your own code:

Contributing your own code
--------------------------

If you are planning to make an extension to or already have a piece of
your own code which could be useful to others, you are very welcome to
contribute it to the community. Before you start making any changes to
the code, you should fork the `espressomd/espresso
<https://github.com/espressomd/espresso>`_ repository and work in a new branch.

It is also generally a good idea to contact the developers on the mailing lists
before you start major coding projects. It might be that someone else is already
working on the problem or has a solution at hand.

You will find more detailed information on our development processes in
`Contributing to ESPResSo
<https://github.com/espressomd/espresso/blob/python/CONTRIBUTING.md>`_.
Please also refer to our `wiki <https://github.com/espressomd/espresso/wiki>`_.

.. _Required Development Tools:

Required Development Tools
^^^^^^^^^^^^^^^^^^^^^^^^^^

-  First of all, please install the dependencies for compiling |es|.
   See the section :ref:`Installation`.

-  To be able to access the development version of |es|, you will need the
   distributed versioning control system git_ and a GitHub account to fork the
   `espressomd/espresso <https://github.com/espressomd/espresso>`_ repository

-  To build the sphinx documentation, you will need the Python packages listed
   in :file:`requirements.txt` in the top-level source directory. To install
   them, use:

   .. code-block:: bash

      pip install --upgrade --user -r requirements.txt

   Note that some distributions now use ``pip`` for Python3 and ``pip2`` for
   Python2.

-  To build the tutorials, you will need LaTeX and Jupyter.

-  To compile the Doxygen code documentation, you will need the tool doxygen_.

All of these tools should be easy to install on most Unix operating systems.

.. _Building the User's guide:

Building the User's guide
-------------------------

If, while reading this documentation, you notice any mistakes or badly (if
at all) described features or commands, you are very welcome to
contribute to the guide and have others benefit from your knowledge.

For this, you should clone the development version at `espressomd/espresso
<https://github.com/espressomd/espresso>`_. Next build the software as shown
in :ref:`Installation` and then the documentation with ``make sphinx``.


.. _git: http://git-scm.com/

.. _doxygen: http://www.doxygen.org/
