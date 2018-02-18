How to get started with SPHINX
==============================

#. Install the latest stable sphinx version via pip::

    pip2 install git+https://github.com/sphinx-doc/sphinx@stable --user --upgrade

#. Install a bibtex extension to Sphinx::

    pip2 install sphinxcontrib-bibtex --user --upgrade

#. Compile the ``sphinx`` target in your build directory (that can take some time
   since we depend on finishing the build of the interface)::

      make sphinx

#. In the directory ``doc/sphinx`` you can find all the source rst files for the user guide.
   You can change those files and rerun ``make sphinx``.

#. When writing a docstring please keep the style to
   `numpy docstrings <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.
