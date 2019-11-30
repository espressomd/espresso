How to get started with SPHINX
==============================

#. Install the latest stable sphinx version and the bibtex extension via pip:

   .. code-block:: bash

       pip3 install --user --upgrade 'sphinx>=1.6.6,!=2.1.0' 'sphinxcontrib-bibtex>=0.3.5'

#. Compile the ``sphinx`` target in your build directory (that can take some time
   since we depend on finishing the build of the interface):

   .. code-block:: bash

         make sphinx

#. In the directory :file:`doc/sphinx` you can find all the source rst files for the user guide.
   You can change those files and rerun ``make sphinx``.

#. When writing a docstring please keep the style to
   `numpy docstrings <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
