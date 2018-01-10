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
#endif

#endif
#endif
