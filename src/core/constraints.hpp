#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include "ParticleRange.hpp"
#include "constraints/Constraint.hpp"
#include "constraints/Constraints.hpp"
#include "energy.hpp"
#include "particle_data.hpp"

namespace Constraints {
extern Constraints<ParticleRange, Constraint> constraints;
}

#endif
