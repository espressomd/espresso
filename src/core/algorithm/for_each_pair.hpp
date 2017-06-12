#ifndef CORE_ALGORITHM_PAIR_LOOP_HPP
#define CORE_ALGORITHM_PAIR_LOOP_HPP

#include <utility>

#include "link_cell.hpp"
#include "verlet_ia.hpp"

namespace Algorithm {
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction, typename VerletCriterion>
void for_each_pair(CellIterator first, CellIterator last,
                   ParticleKernel &&particle_kernel, PairKernel pair_kernel,
                   DistanceFunction &&distance_function,
                   VerletCriterion &&verlet_criterion, bool use_verlet_list,
                   bool rebuild) {
  if (use_verlet_list) {
    verlet_ia(first, last, std::forward<ParticleKernel>(particle_kernel),
              std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function),
              std::forward<VerletCriterion>(verlet_criterion), rebuild);
  } else {
    link_cell(first, last, std::forward<ParticleKernel>(particle_kernel),
              std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function));
  }
}
}

#endif
