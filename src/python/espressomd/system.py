from __future__ import print_function
from espressomd.highlander import highlander

@highlander
class System:
    def __init__(self):
        print("I am a System!")
