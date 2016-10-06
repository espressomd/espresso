from espressomd import script_interface

import espressomd
import sys
import numpy as np


system = espressomd.System()
system.box_l = [10.0, 10.0, 10.0]
system.time_step = 0.01
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.skin = 0.4

for i in range(10):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=23)

h5md = script_interface.PScriptInterface("ScriptInterface::Writer::H5mdScript")
h5md.set_params(filename="test.h5", scriptname=sys.argv[0])
h5md.call_method("init_file")
print h5md.get_params()
h5md.call_method("write")
system.part.add(id=10, pos=np.random.random(3) * system.box_l, type=23)
h5md.call_method("write")
system.part[10].delete()
h5md.call_method("write")
h5md.call_method("close")
