#ifndef INTEGRATORS_VELOCITY_VERLET_NPT_HPP
#define INTEGRATORS_VELOCITY_VERLET_NPT_HPP

#include "config.hpp"

#ifdef NPT
#include "ParticleRange.hpp"
#include "particle_data.hpp"

/** Special propagator for NPT ISOTROPIC
    Propagate the velocities and positions. Integration steps before force
   calculation of the Velocity Verlet integrator: <br> \f[ v(t+0.5 \Delta t) =
   v(t) + 0.5 \Delta t f(t)/m \f] <br> \f[ p(t+\Delta t) = p(t) + \Delta t
   v(t+0.5 \Delta t) \f]

    Propagate pressure, box_length (2 times) and positions, rescale
    positions and velocities and check Verlet list criterion (only NPT) */
void velocity_verlet_npt_step_1(const ParticleRange &particles);

/** Final integration step of the Velocity Verlet+NPT integrator and finalize
    instantaneous pressure calculation:<br>
    \f[ v(t+\Delta t) = v(t+0.5 \Delta t) + 0.5 \Delta t f(t+\Delta t)/m \f] */
void velocity_verlet_npt_step_2(const ParticleRange &particles);
#endif

#endif
