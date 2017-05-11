
#ifndef CORE_CELLITER_H_INCLUDED
#define CORE_CELLITER_H_INCLUDED

#include <utility>
#include "utils/Range.hpp"
#include "ParticleIterator.hpp"
#include "particle_data.hpp"

typedef ParticleList Cell; // Same as in cells.hpp (in order not to introduce
                           // cyclic dependencies)

/** Encapsulates a Cell.
 * This is a helper class to iterate over cells and the particles in them.
 * This class does not take ownership of the encapsulated cell. Warning:
 * Changes to the underlying cells (e.g. local_cells) array might invalidate
 * CellProxies.
 */
struct CellProxy {
  struct NoDomainDecompositionError {};

  /** Constructor.
   * \param c Double pointer to cell
   */
  CellProxy(Cell **c): c(c) {}

  /** Returns a pointer to the cell.
   */
  Cell* cellptr() const { return *c; }
  /** Returns a reference to the cell.
   */
  Cell& cell() const { return **c; }

  using CellParticleIterator = ParticleIterator<Cell **, Particle>;
  /** Returns a begin iterator to the particles of the encapsulated cell.
   */
  CellParticleIterator pbegin() const {
    return CellParticleIterator(c, c + 1, 0);
  }

  /** Returns an end iterator to the particles of the encapsulated cell.
   */
  CellParticleIterator pend() const {
    return CellParticleIterator(c + 1, c + 1, 0);
  }

  /** Returns a range object to iterate over particles of the encapsulates cell.
   * Can be used in range-based for loops.
   */
  Utils::Range<CellParticleIterator> particles() const {
    return Utils::make_range(pbegin(), pend());
  }

  operator Cell*() const { return cellptr(); }
  operator Cell&() const { return cell(); }

  bool operator==(const CellProxy& other) const {
    return cellptr() == other.cellptr();
  }

  bool operator!=(const CellProxy& other) const {
    return !(*this == other);
  }

private:
  Cell **c;
};

/** Iterator over dd.cell_inter.
 * dd.cell_inter[cellidx].nList[i].cell_ind is an index to "cells" and
 * not to local_cells, so it cannot be used to index dd.cell_inter, again.
 * Therefore, this iterator returns only CellProxies and not Indexed ones.
 */
struct CellNeighIterator: public boost::iterator_facade<
                                   CellNeighIterator,
                                   CellProxy,
                                   boost::random_access_traversal_tag,
                                   CellProxy> // Returns value, not reference
{
  /** Constructor
   * \param cellidx Index of the center
   * \poaram i Index of the neighbor
   */
  CellNeighIterator(int cellidx, int i): cellidx(cellidx), i(i) {}

  struct DifferentBaseCellError {};

private:
  const int cellidx;
  int i;

  friend class boost::iterator_core_access;

  value_type dereference() const;

  bool equal(const CellNeighIterator& other) const {
    return cellidx == other.cellidx && i == other.i;
  }

  void increment() { advance(1); }
  void decrement() { advance(-1); }
  void advance(difference_type n) { i += n; }

  difference_type distance_to(const CellNeighIterator& other) const {
    if (cellidx != other.cellidx)
      throw DifferentBaseCellError();
    return i - other.i;
  }
};

/** Encapsulates a Cell.
 * This is a helper class to iterate over cells, the particles in them and
 * their neighborhood. A plain Cell* isn't enough since dd.cell_inter is
 * accessed via a cell index. Therefore, this class not only keeps track of a
 * Cell** but also its index. This class does not take ownership of the
 * encapsulated cell. Warning: Changes to the underlying cells (e.g.
 * local_cells) array might invalidate CellProxies.
 */
struct IndexedCellProxy {
  struct NoDomainDecompositionError {};

  /** Constructor.
   * \param c Double pointer to cell
   * \param idx Index of that cell
   */
  IndexedCellProxy(Cell **c, int idx): c(c), idx(idx) {}

  /** Returns a pointer to the cell.
   */
  Cell* cellptr() const { return *c; }
  /** Returns a reference to the cell.
   */
  Cell& cell() const { return **c; }
  /** Returns the cell index.
   */
  int cellidx() const { return idx; }

  using IdxPtrPair = std::pair<int, Cell*>;
  /** Returns a pair (cell**, index)
   */
  IdxPtrPair idxptrpair() const { return std::make_pair(cellidx(), cellptr()); }

  using CellParticleIterator = ParticleIterator<Cell **, Particle>;
  /** Returns a begin iterator to the particles of the encapsulated cell.
   */
  CellParticleIterator pbegin() const {
    return CellParticleIterator(c, c + 1, 0);
  }

  /** Returns an end iterator to the particles of the encapsulated cell.
   */
  CellParticleIterator pend() const {
    return CellParticleIterator(c + 1, c + 1, 0);
  }

  /** Returns a range object to iterate over particles of the encapsulates cell.
   * Can be used in range-based for loops.
   */
  Utils::Range<CellParticleIterator> particles() const {
    return Utils::make_range(pbegin(), pend());
  }

  /** Returns a range object to iterate over neighbor cells (half-shell).
   * Can be used in range-based for loops.
   * DO NOT use this for ghost cells.
   */
  Utils::Range<CellNeighIterator> hsneigh() const;

  operator Cell*() const { return cellptr(); }
  operator Cell&() const { return cell(); }
  operator int() const { return cellidx(); }
  operator IdxPtrPair() const { return idxptrpair(); }

  bool operator==(const IndexedCellProxy& other) const {
    return cellptr() == other.cellptr();
  }

  bool operator!=(const IndexedCellProxy& other) const {
    return !(*this == other);
  }

private:
  Cell **c;
  int idx;
};


/** Iterator over cells. Do not use this class directly, use local_cells.cells().
 *
 */
struct CellIterator: public boost::iterator_facade<
                              CellIterator,
                              IndexedCellProxy,
                              boost::random_access_traversal_tag,
                              IndexedCellProxy> // Returns value, not reference
{
  typedef IndexedCellProxy value_type;

  /** Constructor
   * \param cells cells array to iterate over (local_cells.cell or ghost_cells.cell)
   * \param i Iteration index
   */
  CellIterator(Cell **cells, int i): cells(cells), i(i) {}

private:
  Cell **cells;
  int i;

  friend class boost::iterator_core_access;

  value_type dereference() const {
    return value_type(&cells[i], i);
  }

  bool equal(const CellIterator& other) const { return i == other.i; }
  void increment() { advance(1); }
  void decrement() { advance(-1); }
  void advance(difference_type n) { i += n; }

  difference_type distance_to(const CellIterator& other) const {
    return i - other.i;
  }
};

inline bool operator==(const IndexedCellProxy& a, const CellProxy& b) {
  return a.cellptr() == b.cellptr();
}

inline bool operator!=(const IndexedCellProxy& a, const CellProxy& b) {
  return !(a == b);
}

#endif
