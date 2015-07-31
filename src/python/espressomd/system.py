from __future__ import print_function
from espressomd.highlander import highlander
from espressomd import _system


@highlander
class System(_system.System):

    def __init__(self):
        _system.System.__init__(self)
