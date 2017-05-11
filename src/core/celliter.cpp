
#include "celliter.hpp"
#include "cells.hpp"
#include "domain_decomposition.hpp"

Utils::Range<CellNeighIterator> IndexedCellProxy::hsneigh() const
{
  if (cell_structure.type != CELL_STRUCTURE_DOMDEC)
    throw NoDomainDecompositionError();
  return Utils::make_range(CellNeighIterator(idx, 0),
                           CellNeighIterator(idx, dd.cell_inter[idx].n_neighbors));
}

CellNeighIterator::value_type CellNeighIterator::dereference() const
{
  return value_type(&dd.cell_inter[cellidx].nList[i].pList);
}

