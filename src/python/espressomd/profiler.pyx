cimport c_profiler


def begin_section(name):
    """
    Start named section in profiler.

    Parameters
    ----------

    name : obj:`str`
        Name of the section
    """

    c_profiler.begin_section(name)


def end_section(name):
    """
    End named section in profiler.

    Parameters
    ----------

    name: obj :`str`
        Name of the section
    """

    c_profiler.end_section(name)
