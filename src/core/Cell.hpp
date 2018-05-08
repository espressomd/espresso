#ifndef CORE_CELL_HPP
#define CORE_CELL_HPP

#include <functional>
#include <vector>

/* Needed for transform_iterator to work with
   lambdas on older compilers. */
#define BOOST_RESULT_OF_USE_DECLTYPE

#include <boost/iterator/transform_iterator.hpp>

#include "particle_data.hpp"
#include "utils/Range.hpp"

class Cell : public ParticleList {
  struct GetReference {
    using result_type = Cell;
    Cell &operator()(std::reference_wrapper<Cell> & cell_ref) const {
      return cell_ref.get();
    }
  };

  using neighbors_t = std::vector<std::reference_wrapper<Cell>>;

public:
  using neighbor_iterator =
      boost::transform_iterator<GetReference, neighbors_t::iterator>;

  /** Topological neighbors of the cell */
  std::vector<std::reference_wrapper<Cell>> m_neighbors;

  /** Interaction pairs */
  std::vector<std::pair<Particle *, Particle *>> m_verlet_list;

  Utils::Range<neighbor_iterator> neighbors() {
    return Utils::make_range(neighbor_iterator(m_neighbors.begin()),
                             neighbor_iterator(m_neighbors.end()));
  }

  void resize(size_t size) {
    realloc_particlelist(static_cast<ParticleList *>(this), this->n = size);
  }

#ifdef LEES_EDWARDS
  int myIndex[3];
#endif
};

#endif
