#ifndef CORE_LB_POPULATION_VIEW_HPP
#define CORE_LB_POPULATION_VIEW_HPP

namespace LB {
namespace MemoryLayout {
/* @brief Tag type indicating that the first index is the mode.
 *
 * E.g. field[i][j] is population i at index j.
 */
struct PopulationIndex;
/* @brief Tag type indicating that the first index is the mode.
 *
 * E.g. field[i][j] is index i and population j.
 */
struct IndexPopulation;
}

template <class T, class Tag = MemoryLayout::PopulationIndex>
struct PopulationView;

template <class T> struct PopulationView<T, MemoryLayout::PopulationIndex> {
private:
  T **const m_field;
  const int m_index;

public:
  explicit PopulationView(T **field, int index)
      : m_field(field), m_index(index) {}

  T &operator()(int pop) { return m_field[pop][m_index]; }
  T const &operator()(int pop) const { return m_field[pop][m_index]; }

  T &operator()(int pop, int offset) { return m_field[pop][m_index + offset]; }
  T const &operator()(int pop, int offset) const { return m_field[pop][m_index + offset]; }
};
}

#endif
