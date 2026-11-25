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
  /** @brief Number of participants residing on this rank. */
  std::uint16_t n_local = 0u;
  /** @brief Participants in creation order, owner first. */
  std::array<int, 4> ids{};

  bool is_free() const { return bond_id < 0; }
  std::span<const int> particle_ids() const { return {ids.data(), n}; }
  BondRecord record() const { return {bond_id, particle_ids()}; }
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
 * participates in serves the per-particle queries.
 */
class BondTable {
public:
  /** @brief Maximal number of particles in a bond. */
  static constexpr std::size_t max_particles = 4u;
  /** @brief Ints per row in the buffers of @ref extract(). */
  static constexpr std::size_t packed_row_size = 2u + max_particles;

private:
  std::vector<BondRow> m_rows;
  /** @brief Free rows, reused by new bonds. */
  std::vector<std::uint32_t> m_free;
  /** @brief For every local particle, the rows it participates in. */
  boost::unordered_flat_map<int, std::vector<std::uint32_t>> m_rows_of;

  inline static std::vector<std::uint32_t> const no_rows{};

  std::vector<std::uint32_t> const &rows_of(int pid) const {
    auto const it = m_rows_of.find(pid);
    return (it == m_rows_of.end()) ? no_rows : it->second;
  }

  /** @brief Find the row of a bond among the rows of a local particle. */
  std::optional<std::uint32_t> find_in(int pid, int bond_id,
                                       std::span<const int> ids) const {
    for (auto const row : rows_of(pid)) {
      if (m_rows[row].record().matches(bond_id, ids)) {
        return row;
      }
    }
    return std::nullopt;
  }

  /**
   * @brief Find the row of a bond among the rows of its local
   * participants, except @p skip_pid.
   */
  std::optional<std::uint32_t> find(int bond_id, std::span<const int> ids,
                                    int skip_pid = -1) const {
    // every local participant lists the bond: search the shortest list
    std::optional<int> best;
    std::size_t best_size = 0u;
    for (auto const pid : ids) {
      if (pid == skip_pid) {
        continue;
      }
      auto const it = m_rows_of.find(pid);
      if (it != m_rows_of.end() and
          (not best or it->second.size() < best_size)) {
        best = pid;
        best_size = it->second.size();
      }
    }
    return best ? find_in(*best, bond_id, ids) : std::nullopt;
  }

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
    r.n_local = 0u;
    r.ids = {};
    std::ranges::copy(ids, r.ids.begin());
    return row;
  }

  void link(int pid, std::uint32_t row) {
    m_rows_of[pid].push_back(row);
    ++m_rows[row].n_local;
  }

  /** @brief Remove a local particle from a row; free the row when no
   *  participant remains on this rank. */
  void unlink(int pid, std::uint32_t row) {
    auto const it = m_rows_of.find(pid);
    if (it == m_rows_of.end()) {
      return;
    }
    auto &rows = it->second;
    auto const pos = std::ranges::find(rows, row);
    if (pos == rows.end()) {
      return; // the particle does not list this row
    }
    rows.erase(pos);
    if (rows.empty()) {
      m_rows_of.erase(it);
    }
    if (--m_rows[row].n_local == 0u) {
      m_rows[row].bond_id = -1;
      m_free.push_back(row);
    }
  }

  /** @brief Remove a row from all its local participants. */
  void drop(std::uint32_t row) {
    auto const ids = m_rows[row].ids;
    auto const n = m_rows[row].n;
    for (std::size_t i = 0u; i < n; ++i) {
      // a particle listed twice in a bond is linked once
      if (std::find(ids.begin(), ids.begin() + i, ids[i]) == ids.begin() + i) {
        unlink(ids[i], row);
      }
    }
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
    for (auto &[pid, rows] : m_rows_of) {
      for (auto &row : rows) {
        row = new_index[row];
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
      link(pid, row);
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
      drop(*row);
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
      drop(row);
    }
    return selected.size();
  }

  /** @brief Forget a particle that no longer resides on this rank. The
   *  bonds of its partners on this rank are kept. */
  void erase(int pid) {
    auto const rows = rows_of(pid);
    for (auto const row : rows) {
      unlink(pid, row);
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
      if (find_in(pid, bond_id, ids)) {
        continue; // the particle had not left
      }
      // the particle is not listed for all its bonds yet: search the
      // bond among the other local participants
      auto const row = find(bond_id, ids, pid);
      link(pid, row ? *row : new_row(bond_id, ids));
    }
    return buffer;
  }
};
