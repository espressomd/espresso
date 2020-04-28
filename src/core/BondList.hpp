/*
 * Copyright (C) 2020 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef ESPRESSO_BONDLIST_HPP
#define ESPRESSO_BONDLIST_HPP

#include <utils/Span.hpp>

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/container/vector.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/version.hpp>
#include <type_traits>

/**
 * @brief Immutable view on a bond.
 *
 * The bond id and ids of the partner particles
 * of a bond can be inspected (but not changed)
 * via this view.
 */
class BondView {
  /* Bond id */
  int m_id = -1;
  Utils::Span<const int> m_partners;

public:
  BondView() = default;
  BondView(int id, Utils::Span<const int> partners)
      : m_id(id), m_partners(partners) {}

  int bond_id() const { return m_id; }
  Utils::Span<const int> const &partner_ids() const { return m_partners; }

  bool operator==(BondView const &rhs) const {
    return m_id == rhs.m_id and boost::equal(m_partners, rhs.m_partners);
  }

  bool operator!=(BondView const &rhs) const { return not(*this == rhs); }
};

/**
 * @brief Bond storage.
 *
 * This holds the bonds of a particle.
 *
 * Implementation notes:
 * Internally the bond list is represented as a sequence of
 * integers. For each bond, first the ids of the bond partners
 * are stored (which are positive numbers), then the negative
 * of the value of the id of the bond plus one (which also is positive) is
 * stored. This way we can use the sign bit of the int as a delimiter,
 * every time a negative number is encountered one bond ends, and
 * a new one starts. This mechanism allows us to efficiently store
 * bonds of different numbers of partners in contiguous memory.
 *
 * This is hidden from the client by providing access to the
 * bonds only thru iterators that restore the semantic meanings
 * of the entries, and know how to proceed to the next bond.
 *
 */
class BondList {
public:
  using storage_type = boost::container::vector<int>;

private:
  using storage_iterator = storage_type::const_iterator;

  /**
   * @brief Find the end of bond.
   *
   * @param it Valid iterator into the underlying storage.
   * @return iterator pointing to the last element of the bond.
   */
  static storage_iterator find_end(storage_iterator it) {
    while (*it >= 0) {
      ++it;
    }
    return it;
  }

  storage_type m_storage;

  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    if (Archive::is_loading::value) {
      size_t size{};
      ar &size;
      m_storage.resize(size);
    }

    if (Archive::is_saving::value) {
      auto size = m_storage.size();
      ar &size;
    }

    ar &boost::serialization::make_array(m_storage.data(), m_storage.size());
  }

public:
  class Iterator
      : public boost::iterator_facade<Iterator, BondView,
                                      boost::forward_traversal_tag, BondView> {
  public:
    explicit Iterator(storage_iterator it) : m_it(it) {}

  private:
    /** Iterator into the bond list */
    storage_iterator m_it;

    friend BondList;
    friend boost::iterator_core_access;
    void increment() { m_it = std::next(find_end(m_it)); }
    bool equal(Iterator const &other) const { return this->m_it == other.m_it; }
    BondView dereference() const {
      auto const id_pos = find_end(m_it);
      auto const partners_begin = m_it;
      auto const partners_end = id_pos;
      return {-(*id_pos) - 1,
              Utils::make_span(std::addressof(*partners_begin),
                               std::distance(partners_begin, partners_end))};
    }
  };

public:
  using value_type = BondView;
  using reference = std::add_lvalue_reference_t<BondView>;
  using const_reference = std::add_const_t<reference>;
  using size_type = size_t;
  using difference_type = ptrdiff_t;
  using iterator = Iterator;
  using const_iterator = Iterator;

  BondList() = default;
  BondList(BondList const &) = default;
  BondList(BondList &&) = default;
  BondList &operator=(BondList const &rhs) {
    if (this != std::addressof(rhs)) {
      m_storage = rhs.m_storage;
    }

    return *this;
  }
  // NOLINTNEXTLINE(bugprone-exception-escape)
  BondList &operator=(BondList &&rhs) noexcept {
    if (this != std::addressof(rhs)) {
      std::swap(m_storage, rhs.m_storage);
    }

    return *this;
  }

  /**
   * @brief Iterator to the beginning of the range of bonds in the list.
   */
  const_iterator begin() const { return Iterator{m_storage.begin()}; }
  /**
   * @brief Iterator past the end of the range of bonds in the list.
   */
  const_iterator end() const { return Iterator{m_storage.end()}; }

  /**
   * @brief Add a bond to the list.
   *
   * @param bond Bond to add.
   */
  void insert(BondView const &bond) {
    boost::copy(bond.partner_ids(), std::back_inserter(m_storage));
    assert(bond.bond_id() >= 0);
    m_storage.push_back(-(bond.bond_id() + 1));
  }

  /**
   * @brief Erase a bond from the list.
   * @param pos Iterator pointing to element to erase.
   * @return iterator pointing one past the erased element.
   */
  const_iterator erase(const_iterator pos) {
    return Iterator{m_storage.erase(pos.m_it, std::next(find_end(pos.m_it)))};
  }

  /**
   * @brief Number of bonds.
   * @return The number of bonds in the list.
   */
  size_type size() const { return std::distance(begin(), end()); }

  /**
   * @brief Erase all bonds from the list.
   */
  void clear() { m_storage.clear(); }

  /**
   * @brief Check if the are any bonds in the list.
   */
  bool empty() const { return m_storage.empty(); }

  // NOLINTNEXTLINE(bugprone-exception-escape)
  friend void swap(BondList &lhs, BondList &rhs) {
    using std::swap;

    swap(lhs.m_storage, rhs.m_storage);
  }
};

/**
 * @brief Check if there is a specific bond in a bond list.
 *
 * @param bonds List of bonds to check
 * @param partner_id Id of the bond partner
 * @param bond_id Id of the bond parameters.
 * @return True iff there is a bond to the specified id of the specified type.
 */
inline bool pair_bond_exists_on(BondList const &bonds, int partner_id,
                                int bond_id) {
  return boost::algorithm::any_of(bonds, [=](BondView const &bond) {
    return (bond.bond_id() == bond_id) and
           (bond.partner_ids()[0] == partner_id);
  });
}

#endif // ESPRESSO_BONDLIST_HPP
