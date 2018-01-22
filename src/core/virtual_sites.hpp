#ifndef VIRTUAL_SITES_HPP
#define VIRTUAL_SITES_HPP

#include "config.hpp"

#ifdef VIRTUAL_SITES
#include "virtual_sites/VirtualSites.hpp" 

/** @brief get active virtual sites implementation */
const std::shared_ptr<VirtualSites>& virtual_sites();

/** @brief Set active virtual sites implementation */
void set_virtual_sites(std::shared_ptr<VirtualSites> const& v);

#ifdef VIRTUAL_SITES_RELATIVE
int vs_relate_to(int part_num, int relate_to);

// Setup the virtual_sites_relative properties of a particle so that the given virtaul particle will follow the given real particle
// Local version, expects both particles to be accessible through local_particles
// and only executes the changes on the virtual site locally
int local_vs_relate_to(int part_num, int relate_to);


#endif
#endif
#endif
