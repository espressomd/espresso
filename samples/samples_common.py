#!/usr/bin/env python

"""
Module for methods common to several sample files.
"""

import tempfile


def open(*args, **kwargs):
    """
    Return a temporary file object.

    Parameters
    ----------

    *args: anytype
           All arguments are discarded

    **args: anytype
            All keyword arguments are discarded

    Returns
    -------
    file-like object
    """
    return tempfile.TemporaryFile()
