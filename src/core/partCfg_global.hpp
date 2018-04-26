#ifndef CORE_PART_CFG_GLOBAL_HPP
#define CORE_PART_CFG_GLOBAL_HPP
 
#include "PartCfg.hpp"

/**
 * @brief Particles' current configuration.
 *
 * Particle coordinates are unfolded.
 * For documentation see @class ParticleCache
 */
PartCfg & partCfg(std::unique_ptr<PartCfg> init = std::unique_ptr<PartCfg>{});

#endif
