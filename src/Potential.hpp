/*
 * Force.hpp
 *
 *  Created on: May 30, 2014
 *      Author: olenz
 */

#ifndef _POTENTIAL_HPP
#define _POTENTIAL_HPP

#include "SystemInterface.hpp"

/**
 * Generic abstract potential class.
 * All potentials should be derived from this one.
 */
class Potential {
public:
	virtual void computeForces(SystemInterface &s) = 0;
};

#endif /* _POTENTIAL_HPP */
