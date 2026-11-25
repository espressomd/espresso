/*
 * Copyright (C) 2026 The ESPResSo project
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

#pragma once

#include <boost/unordered/unordered_flat_map.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

/**
 * @brief Read-only view on one bond.
 *
 * A bond is the key of a bonded interaction and the ids of the
 * participating particles in creation order. The first particle is the
 * owner of the bond: the bond is evaluated on the MPI rank where the owner
 * resides.
 */
class BondRecord {
  int m_bond_id = -1;
  std::span<const int> m_ids;

public:
  BondRecord() = default;
  BondRecord(int bond_id, std::span<const int> particle_ids)
      : m_bond_id(bond_id), m_ids(particle_ids) {}

  /** @brief Key of the bonded interaction. */
  int bond_id() const { return m_bond_id; }
  /** @brief Participants in creation order; the first one is the owner. */
  std::span<const int> particle_ids() const { return m_ids; }
  int owner_id() const { return m_ids.front(); }
  std::span<const int> partner_ids() const { return m_ids.subspan(1u); }
  bool has_participant(int pid) const {
    return std::ranges::find(m_ids, pid) != m_ids.end();
  }
  /** @brief Same bond id and same participants in the same order. */
  bool matches(int bond_id, std::span<const int> particle_ids) const {
    return m_bond_id == bond_id and std::ranges::equal(m_ids, particle_ids);
  }
};

/** @brief One bond in a @ref BondTable. */
struct BondRow {
  /** @brief Key of the bonded interaction, negative for a free row. */
  int bond_id = -1;
  /** @brief Number of participants. */
  std::uint16_t n = 0u;
  /** @brief Bit @c i is set if participant @c i resides on this rank. */
  std::uint16_t local_mask = 0u;
  /** @brief Participants in creation order, owner first. */
  std::array<int, 4> ids{};

  bool is_free() const { return bond_id < 0; }
  std::span<const int> particle_ids() const { return {ids.data(), n}; }
  BondRecord record() const { return {bond_id, particle_ids()}; }
  /** @brief Bits of the positions at which a particle participates. */
  std::uint16_t mask_of(int pid) const {
    std::uint16_t mask = 0u;
    for (std::size_t i = 0u; i < n; ++i) {
      if (ids[i] == pid) {
        mask = static_cast<std::uint16_t>(mask | (1u << i));
      }
    }
    return mask;
  }
  /** @brief Whether a participant resides on this rank. */
  bool is_local(int pid) const { return (local_mask & mask_of(pid)) != 0u; }
};

/**
 * @brief The bonds of the particles residing on this MPI rank.
 *
 * Every bond with at least one participant on this rank is stored once, as
 * a fixed-size row. A bond is therefore known on every rank where one of
 * its participants resides, and on no other rank. The bond information is
 * kept separate from the particles: when a particle moves to another rank,
 * the rows it participates in are shipped with it (@ref extract(),
 * @ref insert()).
 *
 * A bond is identified by its bond id and its participants in creation
 * order, and is stored at most once: e.g. (1, 2) and (2, 1) are different
 * bonds, but (1, 2) cannot be stored twice.
 *
 * The rows are scanned linearly to build the flat bond lists of the force
 * kernels (@ref rows()). An index from every local particle to the rows it
 * participates in serves the per-particle queries. A second index lists,
 * for every participant residing on another rank, the rows it appears in,
 * so that all bonds of a particle can be removed on every rank without
 * communication (@ref remove_bonds_of()). Only bonds across rank
 * boundaries have entries in that second index.
 */
class BondTable {
public:
  /** @brief Maximal number of particles in a bond. */
  static constexpr std::size_t max_particles = 4u;
  /** @brief Ints per row in the buffers of @ref extract(). */
  static constexpr std::size_t packed_row_size = 2u + max_particles;

private:
  using Index = boost::unordered_flat_map<int, std::vector<std::uint32_t>>;

  std::vector<BondRow> m_rows;
  /** @brief Free rows, reused by new bonds. */
  std::vector<std::uint32_t> m_free;
  /** @brief For every local particle, the rows it participates in, in the
   *  order it became a local participant of them. */
  Index m_rows_of;
  /** @brief For every particle residing on another rank, the rows it
   *  participates in. */
  Index m_remote_rows_of;

  inline static std::vector<std::uint32_t> const no_rows{};

  std::vector<std::uint32_t> const &rows_of(int pid) const {
    auto const it = m_rows_of.find(pid);
    return (it == m_rows_of.end()) ? no_rows : it->second;
  }

  /** @brief Call @p f once for every distinct participant of a row. */
  template <typename F> void for_each_distinct(BondRow const &r, F &&f) const {
    for (std::size_t i = 0u; i < r.n; ++i) {
      if (std::find(r.ids.begin(), r.ids.begin() + i, r.ids[i]) ==
          r.ids.begin() + i) {
        f(r.ids[i]);
      }
    }
  }

  static void push(Index &index, int pid, std::uint32_t row) {
    index[pid].push_back(row);
  }

  static void pop(Index &index, int pid, std::uint32_t row) {
    auto const it = index.find(pid);
    if (it == index.end()) {
      return;
    }
    auto &rows = it->second;
    // recently added rows are the most likely to be removed
    auto const pos = std::find(rows.rbegin(), rows.rend(), row);
    if (pos == rows.rend()) {
      return;
    }
    rows.erase(std::next(pos).base());
    if (rows.empty()) {
      index.erase(it);
    }
  }

  /** @brief Find the row of a bond among the rows of its participants. */
  std::optional<std::uint32_t> find(int bond_id,
                                    std::span<const int> ids) const {
    // every participant of a stored bond lists it in exactly one of the two
    // indices: search the participant with the fewest rows
    auto const &remote_rows_of = [this](int pid) -> auto const & {
      auto const it = m_remote_rows_of.find(pid);
      return (it == m_remote_rows_of.end()) ? no_rows : it->second;
    };
    if (ids.empty()) {
      return std::nullopt;
    }
    auto best = ids[0];
    auto best_size = std::numeric_limits<std::size_t>::max();
    for (auto const pid : ids) {
      auto const size = rows_of(pid).size() + remote_rows_of(pid).size();
      if (size < best_size) {
        best = pid;
        best_size = size;
      }
    }
    for (auto const *rows : {&rows_of(best), &remote_rows_of(best)}) {
      for (auto const row : *rows) {
        if (m_rows[row].record().matches(bond_id, ids)) {
          return row;
        }
      }
    }
    return std::nullopt;
  }

  /** @brief Store a bond whose participants all reside on other ranks;
   *  see @ref set_local(). */
  std::uint32_t new_row(int bond_id, std::span<const int> ids) {
    assert(bond_id >= 0);
    assert(not ids.empty() and ids.size() <= max_particles);
    std::uint32_t row;
    if (m_free.empty()) {
      row = static_cast<std::uint32_t>(m_rows.size());
      m_rows.emplace_back();
    } else {
      row = m_free.back();
      m_free.pop_back();
    }
    auto &r = m_rows[row];
    r.bond_id = bond_id;
    r.n = static_cast<std::uint16_t>(ids.size());
    r.local_mask = 0u;
    r.ids = {};
    std::ranges::copy(ids, r.ids.begin());
    for_each_distinct(r, [&](int pid) { push(m_remote_rows_of, pid, row); });
    return row;
  }

  /** @brief Mark a participant of a row as residing on this rank. */
  void set_local(int pid, std::uint32_t row) {
    auto &r = m_rows[row];
    auto const mask = r.mask_of(pid);
    if (mask == 0u or (r.local_mask & mask) != 0u) {
      return; // not a participant, or already local
    }
    r.local_mask = static_cast<std::uint16_t>(r.local_mask | mask);
    pop(m_remote_rows_of, pid, row);
    push(m_rows_of, pid, row);
  }

  /** @brief Mark a participant of a row as residing on another rank; free
   *  the row when no participant remains on this rank. */
  void clear_local(int pid, std::uint32_t row) {
    auto &r = m_rows[row];
    auto const mask = r.mask_of(pid);
    if ((r.local_mask & mask) == 0u) {
      return; // the particle is not a local participant of this row
    }
    r.local_mask = static_cast<std::uint16_t>(r.local_mask & ~mask);
    pop(m_rows_of, pid, row);
    if (r.local_mask == 0u) {
      free_row(row);
    } else {
      push(m_remote_rows_of, pid, row);
    }
  }

  /** @brief Remove a row from the indices of all its participants. */
  void free_row(std::uint32_t row) {
    auto &r = m_rows[row];
    for_each_distinct(r, [&](int pid) {
      pop(r.is_local(pid) ? m_rows_of : m_remote_rows_of, pid, row);
    });
    r.bond_id = -1;
    r.local_mask = 0u;
    m_free.push_back(row);
  }

public:
  /** @brief All bonds a local particle participates in, owned or not, in
   *  insertion order. */
  auto view(int pid) const {
    return rows_of(pid) | std::views::transform([this](std::uint32_t row) {
             return m_rows[row].record();
           });
  }

  /** @brief The bonds a local particle owns, in insertion order. */
  auto owned(int pid) const {
    return view(pid) | std::views::filter([pid](BondRecord const &bond) {
             return bond.owner_id() == pid;
           });
  }

  /** @brief Whether the bond is stored. */
  bool contains(int bond_id, std::span<const int> ids) const {
    return find(bond_id, ids).has_value();
  }

  /** @brief All rows, including free ones (see @ref BondRow::is_free()). */
  std::span<const BondRow> rows() const { return m_rows; }

  /** @brief Number of bonds. */
  std::size_t size() const { return m_rows.size() - m_free.size(); }
  bool empty() const { return size() == 0u; }
  /** @brief Number of local particles with bonds. */
  std::size_t n_particles() const { return m_rows_of.size(); }

  void clear() {
    m_rows.clear();
    m_free.clear();
    m_rows_of.clear();
    m_remote_rows_of.clear();
  }

  /** @brief Whether most rows are free, see @ref compact(). */
  bool is_wasteful() const {
    return m_free.size() > 1024u and 2u * m_free.size() > m_rows.size();
  }

  /** @brief Drop the free rows, which keeps the order of the bonds. */
  void compact() {
    std::vector<std::uint32_t> new_index(m_rows.size());
    std::uint32_t n_live = 0u;
    for (std::size_t row = 0u; row < m_rows.size(); ++row) {
      if (not m_rows[row].is_free()) {
        new_index[row] = n_live;
        m_rows[n_live++] = m_rows[row];
      }
    }
    m_rows.resize(n_live);
    m_rows.shrink_to_fit();
    m_free.clear();
    for (auto *index : {&m_rows_of, &m_remote_rows_of}) {
      for (auto &[pid, rows] : *index) {
        for (auto &row : rows) {
          row = new_index[row];
        }
      }
    }
  }

  /**
   * @brief Store a bond.
   *
   * Nothing is stored if the bond is stored already.
   *
   * @param bond_id     Key of the bonded interaction.
   * @param ids         Participants in creation order, owner first.
   * @param local_pids  The distinct participants residing on this rank.
   */
  void add(int bond_id, std::span<const int> ids,
           std::span<const int> local_pids) {
    if (local_pids.empty() or contains(bond_id, ids)) {
      return;
    }
    auto const row = new_row(bond_id, ids);
    for (auto const pid : local_pids) {
      set_local(pid, row);
    }
  }

  /**
   * @brief Remove a bond, given by its bond id and its participants in
   * creation order.
   * @return whether the bond was stored.
   */
  bool remove(int bond_id, std::span<const int> ids) {
    auto const row = find(bond_id, ids);
    if (row) {
      free_row(*row);
    }
    return row.has_value();
  }

  /**
   * @brief Remove the bonds of a local particle accepted by a predicate.
   * @param pid   Local particle.
   * @param pred  Callable <tt>bool(BondRecord const &)</tt>.
   * @return number of removed bonds.
   */
  template <typename Pred> std::size_t remove_if(int pid, Pred &&pred) {
    std::vector<std::uint32_t> selected;
    for (auto const row : rows_of(pid)) {
      if (pred(m_rows[row].record())) {
        selected.push_back(row);
      }
    }
    for (auto const row : selected) {
      free_row(row);
    }
    return selected.size();
  }

  /**
   * @brief Remove every bond a particle participates in, whether it resides
   * on this rank or not. Calling this on all ranks removes the bonds of the
   * particle everywhere, without communication.
   * @return number of removed bonds.
   */
  std::size_t remove_bonds_of(int pid) {
    std::vector<std::uint32_t> rows = rows_of(pid);
    if (auto const it = m_remote_rows_of.find(pid);
        it != m_remote_rows_of.end()) {
      rows.insert(rows.end(), it->second.begin(), it->second.end());
    }
    for (auto const row : rows) {
      free_row(row);
    }
    return rows.size();
  }

  /** @brief Forget a particle that no longer resides on this rank. The
   *  bonds of its partners on this rank are kept. */
  void erase(int pid) {
    auto const rows = rows_of(pid);
    for (auto const row : rows) {
      clear_local(pid, row);
    }
  }

  /**
   * @brief Move the bonds of a particle into a buffer, e.g. to send them
   * along with the particle to another rank, and forget the particle.
   * Appends the number of bonds, then per bond
   * <tt>[bond_id, n, id_0, ..., id_3]</tt>.
   */
  void extract(int pid, std::vector<int> &buffer) {
    copy_to(pid, buffer);
    erase(pid);
  }

  /** @brief Append the bonds of a particle to a buffer, in the format of
   *  @ref extract(). */
  void copy_to(int pid, std::vector<int> &buffer) const {
    auto const &rows = rows_of(pid);
    buffer.push_back(static_cast<int>(rows.size()));
    for (auto const row : rows) {
      auto const &r = m_rows[row];
      buffer.push_back(r.bond_id);
      buffer.push_back(r.n);
      buffer.insert(buffer.end(), r.ids.begin(), r.ids.end());
    }
  }

  /**
   * @brief Store the bonds of a particle that arrived on this rank, read
   * from a buffer filled by @ref extract(). Bonds already stored, because
   * a partner resides on this rank, are shared.
   * @return the remainder of the buffer.
   */
  std::span<const int> insert(int pid, std::span<const int> buffer) {
    auto const malformed = [] {
      throw std::runtime_error("Malformed bond data");
    };
    if (buffer.empty() or buffer[0] < 0) {
      malformed();
    }
    auto const n_bonds = static_cast<std::size_t>(buffer[0]);
    buffer = buffer.subspan(1u);
    if (buffer.size() / packed_row_size < n_bonds) {
      malformed();
    }
    for (std::size_t i = 0u; i < n_bonds; ++i) {
      auto const bond_id = buffer[0];
      if (bond_id < 0 or buffer[1] < 1 or
          static_cast<std::size_t>(buffer[1]) > max_particles) {
        malformed();
      }
      auto const ids = buffer.subspan(2u, static_cast<std::size_t>(buffer[1]));
      buffer = buffer.subspan(packed_row_size);
      // the bond is stored already if the particle had not left, or if a
      // partner resides on this rank
      auto const row = find(bond_id, ids);
      set_local(pid, row ? *row : new_row(bond_id, ids));
    }
    return buffer;
  }
};
