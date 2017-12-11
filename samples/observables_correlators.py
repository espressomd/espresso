# This scripts demonstrates the measurement of the mean square displacement
# using the Observables/Correlators mechanism

from espressomd import System
from espressomd.observables import *
from espressomd.correlators import *
from numpy import *


# System setup
s = System()
s.part.add(pos=(0, 0, 0), v=(1, 2, 3))
s.time_step = 0.01
s.cell_system.skin = 0
s.cell_system.set_n_square(use_verlet_lists=False)
s.thermostat.set_langevin(kT=1, gamma=10)
s.integrator.run(1000)

# Initialize obzervable for a particle with id 0
p = ParticlePositions(ids=(0,))
# ASk the observable for its parameters
print p.get_params()
# Calculate and return current value
print p.calculate()
# Return stored current value
print p.value()


# Instance a correlator correlating the p observable with itself, calculating the mean squared displacement (msd).
c = Correlator(tau_lin=16, tau_max=1000, dt=0.01, obs1=p,
               corr_operation="square_distance_componentwise", compress1="discard1")
# Instance a correlator calculating the FCS autocorrelation function from particle positions, using the symmetric focal spot with wx=wy=wz=10 (sigma)
fcs = Correlator(tau_lin=16, tau_max=10000, dt=0.1, obs1=p,
                 corr_operation="fcs_acf", args=[10, 10, 10], compress1="discard2")
# Ask the correlator for its parameters
print c.get_params()

# Register the correlator for auto updating at the interval given by its dt (currently every timestep)
s.auto_update_correlators.add(c)
s.auto_update_correlators.add(fcs)

# Integrate
s.integrator.run(300000)

# Finalize the correlation calculation and write the results to a file
c.finalize()
savetxt("res.dat", c.result())
fcs.finalize()
savetxt("fcs.dat", fcs.result())
