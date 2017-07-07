#ifndef CORE_PARTICLE_ITERATOR_HPP
#define CORE_PARTICLE_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>

template <typename BidirectionalIterator, typename Particle>
struct ParticleIterator : public boost::iterator_facade<
                              ParticleIterator<BidirectionalIterator, Particle>,
                              Particle, boost::forward_traversal_tag> {
  ParticleIterator(BidirectionalIterator cell, BidirectionalIterator end,
                   int part_id)
      : m_cell(cell), m_end(end), m_part_id(part_id) {
    /* Jump to first actual particle */
    if ((m_cell != m_end) && (*m_cell)->n == 0) {
      increment();
    }
  }

private:
  friend class boost::iterator_core_access;

  void increment() {
    /* If we are not at the end of the cells,
       there actually are particles in this cells
       and if we are not at the last particle in the
       cell we can just increment the particle id.
    */
    if ((m_cell != m_end) && ((*m_cell)->n > 0) &&
        (m_part_id < ((*m_cell)->n - 1))) {
      /* Next part in same cell */
      ++m_part_id;
    } else {
      m_part_id = 0;

      /* Don't run over the end */
      if (m_cell != m_end)
        ++m_cell;

      /* Find next cell with particles */
      while ((m_cell != m_end) && ((*m_cell)->n == 0)) {
        ++m_cell;
      }
    }
  }

  bool equal(ParticleIterator const &rhs) const {
    return (m_cell == (rhs.m_cell)) && (m_part_id == rhs.m_part_id);
  }

  Particle &dereference() const { return (*m_cell)->part[m_part_id]; }

  BidirectionalIterator m_cell, m_end;
  int m_part_id;
};

#endif
