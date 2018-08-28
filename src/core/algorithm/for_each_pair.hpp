#ifndef CORE_ALGORITHM_PAIR_LOOP_HPP
#define CORE_ALGORITHM_PAIR_LOOP_HPP

#include <utility>

#include "link_cell.hpp"
#include "verlet_ia.hpp"

namespace Algorithm {
/**
 * @brief Run single and pair kernel for each particle (pair) from cell range.
 *
 * Iterates over all cells in [first, last), and calls particle_kernel for
 * each particle in the cells. Then, for every particle pair within the
 * cell and for each pair with the cells neighbors, distance_function is
 * evaluated and verlet_criterion is evaluated with the calculated distance.
 * Iff true, the pair_kernel is called.
 *
 * For details see verlet_ia and link_cell.
 *
 * Requirements on the types:
 * The Cell type has to provide a function neighbors() that returns
 * a cell range coprised of the topological neighbors of the cell,
 * excluding the cell itself. The cells have to provide a m_verlet_list
 * container that can be used to store particle pairs. It can be empty and is
 * not touched if use_verlet_list is false.
 *
 * verlet_criterion(p1, p2, distance_function(p1, p2)) has to be valid and
 * convertible to bool.
 *
 * ParticleKernel has to provide an ::operator() that can be called with a
 * particle reference.
 * PairKernel has to provide an ::operator() that can be called with two
 * particle references an a distance.
 */
template <typename CellIterator, typename ParticleKernel, typename PairKernel,
          typename DistanceFunction, typename VerletCriterion>
void for_each_pair(CellIterator first, CellIterator last,
                   ParticleKernel &&particle_kernel, PairKernel &&pair_kernel,
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
} // namespace Algorithm

#endif
