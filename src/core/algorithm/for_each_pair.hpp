#ifndef CORE_ALGORITHM_PAIR_LOOP_HPP
#define CORE_ALGORITHM_PAIR_LOOP_HPP

#include <utility>

#include "cells.hpp"
#include "link_cell.hpp"
#include "utils/Batch.hpp"
#include "verlet_ia.hpp"

namespace Utils {
struct True {
  template <typename... T> bool operator()(T...) const { return true; }
};
}

namespace Algorithm {
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction, typename VerletCriterion = Utils::True>
void for_each_pair(CellIterator first, CellIterator last,
                   ParticleKernel &&particle_kernel, PairKernel pair_kernel,
                   DistanceFunction &&distance_function,
                   VerletCriterion &&verlet_criterion = VerletCriterion{}) {
  if (cell_structure.use_verlet_list) {
    verlet_ia(first, last, std::forward<ParticleKernel>(particle_kernel),
              std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function),
              std::forward<VerletCriterion>(verlet_criterion));
  } else {
    link_cell(first, last, std::forward<ParticleKernel>(particle_kernel),
              std::forward<PairKernel>(pair_kernel),
              std::forward<DistanceFunction>(distance_function));
  }
}
}

#endif
