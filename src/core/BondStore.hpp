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
#include <span>
#include <vector>

/**
 * @brief Read-only view on one bond record.
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

/**
 * @brief Read-only range over bond records stored contiguously as
 * <tt>[bond_id, n, id_0, ..., id_{n-1}]</tt>.
 */
class BondRecordRange {
  std::span<const int> m_data;

public:
  class const_iterator {
    int const *m_ptr = nullptr;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = BondRecord;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = BondRecord;

    const_iterator() = default;
    explicit const_iterator(int const *ptr) : m_ptr(ptr) {}

    BondRecord operator*() const {
      return {m_ptr[0], std::span<const int>(
                            m_ptr + 2, static_cast<std::size_t>(m_ptr[1]))};
    }
    const_iterator &operator++() {
      m_ptr += 2 + m_ptr[1];
      return *this;
    }
    const_iterator operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }
    bool operator==(const_iterator const &rhs) const {
      return m_ptr == rhs.m_ptr;
    }
  };

  BondRecordRange() = default;
  explicit BondRecordRange(std::span<const int> data) : m_data(data) {}

  const_iterator begin() const { return const_iterator{m_data.data()}; }
  const_iterator end() const {
    return const_iterator{m_data.data() + m_data.size()};
  }
  bool empty() const { return m_data.empty(); }
  /** @brief Number of records (linear time). */
  std::size_t size() const {
    return static_cast<std::size_t>(std::distance(begin(), end()));
  }
  /** @brief The stored ints, as laid out in memory. */
  std::span<const int> data() const { return m_data; }
};

/**
 * @brief Bonds of the particles residing on this MPI rank.
 *
 * Every particle participates in some bonds; for each local particle, the
 * store holds a full copy of each of these bonds, whether the particle
 * owns them or not. A bond is therefore known on every MPI rank where one
 * of its participants resides, and on no other rank. The bond information
 * is kept separate from the particles: when a particle moves to another
 * rank, its records are shipped with it (see @ref extract() and
 * @ref insert()).
 *
 * The records of all particles are stored in one contiguous buffer, one
 * block per particle, as <tt>[bond_id, n, id_0, ..., id_{n-1}]</tt>. Within
 * a block, the records of the bonds the particle owns come first. Changes
 * relocate a block to the end of the buffer when it grows; @ref compact()
 * rewrites the buffer in the order of the local particles, which keeps the
 * bonds of neighboring particles next to each other in memory and gives
 * access by the index of the particle in that order, without lookups.
 *
 * Read-only functions are safe to call concurrently.
 */
class BondStore {
public:
  /** @brief Maximal number of particles in a bond. */
  static constexpr std::size_t max_particles = 4u;

  /** @brief Number of owned bonds, by number of partners. */
  struct Counts {
    int pair = 0;
    int angle = 0;
    int dihedral = 0;
  };

private:
  static constexpr std::size_t header_size = 2u;

  struct Block {
    /** @brief Position of the block in @c m_data. */
    std::uint32_t offset = 0u;
    /** @brief Number of ints in use. */
    std::uint32_t size = 0u;
    /** @brief Number of ints reserved for the block. */
    std::uint32_t capacity = 0u;
    /** @brief Number of ints of the owned records, at the front. */
    std::uint32_t owned_size = 0u;
    Counts counts{};
  };

  std::vector<int> m_data;
  boost::unordered_flat_map<int, Block> m_blocks;
  /** @brief Ints reserved in @c m_data but no longer used by any block. */
  std::size_t m_garbage = 0u;
  /** @brief Ints in use by the blocks. */
  std::size_t m_live = 0u;
  /** @brief Blocks in the order of the last @ref compact(). */
  std::vector<Block> m_ordered;
  bool m_ordered_valid = false;

  static void update_counts(Counts &counts, std::size_t n_particles,
                            int delta) {
    switch (n_particles) {
    case 2u:
      counts.pair += delta;
      break;
    case 3u:
      counts.angle += delta;
      break;
    case 4u:
      counts.dihedral += delta;
      break;
    default:
      break;
    }
  }

  static std::uint32_t record_size(int const *record) {
    return static_cast<std::uint32_t>(header_size) +
           static_cast<std::uint32_t>(record[1]);
  }

  Block const *find(int pid) const {
    auto const it = m_blocks.find(pid);
    return (it == m_blocks.end()) ? nullptr : &it->second;
  }

  std::span<const int> span_of(Block const &block, std::uint32_t size) const {
    return {m_data.data() + block.offset, size};
  }

  /** @brief Make room for @p n more ints in the block, at its end. */
  Block &grow(int pid, std::uint32_t n) {
    if (m_garbage > 1024u and 2u * m_garbage > m_data.size()) {
      compact({}); // reclaim unused ints, keeping the current order
    }
    auto &block = m_blocks[pid];
    if (block.size + n <= block.capacity) {
      return block;
    }
    auto const capacity = std::max(block.size + n, 2u * block.capacity);
    auto const offset = static_cast<std::uint32_t>(m_data.size());
    m_data.resize(m_data.size() + capacity);
    std::copy_n(m_data.begin() + block.offset, block.size,
                m_data.begin() + offset);
    m_garbage += block.capacity;
    block.offset = offset;
    block.capacity = capacity;
    return block;
  }

  /** @brief Remove @p length ints at @p pos of a block, shifting the rest. */
  void erase_in_block(Block &block, std::uint32_t pos, std::uint32_t length) {
    auto const first = m_data.begin() + block.offset;
    std::copy(first + pos + length, first + block.size, first + pos);
    block.size -= length;
    m_live -= length;
    if (pos < block.owned_size) {
      block.owned_size -= length;
    }
  }

  void changed() { m_ordered_valid = false; }

public:
  /** @brief All records of a particle, owned ones first. */
  BondRecordRange view(int pid) const {
    auto const *block = find(pid);
    return block ? BondRecordRange(span_of(*block, block->size))
                 : BondRecordRange{};
  }

  /** @brief The records of the bonds a particle owns, in insertion order. */
  BondRecordRange owned(int pid) const {
    auto const *block = find(pid);
    return block ? BondRecordRange(span_of(*block, block->owned_size))
                 : BondRecordRange{};
  }

  /** @brief Number of bonds a particle owns, in constant time. */
  Counts owned_counts(int pid) const {
    auto const *block = find(pid);
    return block ? block->counts : Counts{};
  }

  /** @brief Whether the particle stores a record with exactly these
   *  participants. */
  bool contains(int pid, int bond_id, std::span<const int> particle_ids) const {
    return std::ranges::any_of(view(pid), [&](BondRecord const &bond) {
      return bond.matches(bond_id, particle_ids);
    });
  }

  /** @brief Whether any particle stores records. */
  bool empty() const { return m_blocks.empty(); }
  /** @brief Number of particles with records. */
  std::size_t n_particles() const { return m_blocks.size(); }
  /** @brief Number of ints allocated, including unused ones. */
  std::size_t capacity() const { return m_data.size(); }
  /** @brief Number of ints in use. */
  std::size_t size_in_use() const { return m_live; }
  /** @brief Whether @ref compact() would free a significant amount of
   *  memory. */
  bool is_wasteful() const {
    return m_data.size() > m_live + m_live / 4u + 4096u;
  }

  void clear() {
    m_data.clear();
    m_blocks.clear();
    m_ordered.clear();
    m_garbage = 0u;
    m_live = 0u;
    changed();
  }

  /**
   * @brief Store a bond record in a particle.
   * @param pid           Particle that stores the record.
   * @param bond_id       Key of the bonded interaction.
   * @param particle_ids  Participants in creation order, owner first.
   */
  void add(int pid, int bond_id, std::span<const int> particle_ids) {
    assert(bond_id >= 0);
    assert(not particle_ids.empty() and particle_ids.size() <= max_particles);
    auto const length =
        static_cast<std::uint32_t>(header_size + particle_ids.size());
    auto &block = grow(pid, length);
    auto const owned = particle_ids.front() == pid;
    auto const first = m_data.begin() + block.offset;
    auto const pos = owned ? block.owned_size : block.size;
    // shift the records that follow the insertion point
    std::copy_backward(first + pos, first + block.size,
                       first + block.size + length);
    first[pos] = bond_id;
    first[pos + 1] = static_cast<int>(particle_ids.size());
    std::ranges::copy(particle_ids, first + pos + header_size);
    block.size += length;
    m_live += length;
    if (owned) {
      block.owned_size += length;
      update_counts(block.counts, particle_ids.size(), +1);
    }
    changed();
  }

  /**
   * @brief Remove from a particle the first record with exactly the given
   * bond id and participants, in the given order.
   * @return whether a record was removed.
   */
  bool remove(int pid, int bond_id, std::span<const int> particle_ids) {
    auto const it = m_blocks.find(pid);
    if (it == m_blocks.end()) {
      return false;
    }
    auto &block = it->second;
    auto const *data = m_data.data() + block.offset;
    for (std::uint32_t pos = 0u; pos < block.size;) {
      auto const size = record_size(data + pos);
      if (BondRecord(data[pos], {data + pos + header_size, size - 2u})
              .matches(bond_id, particle_ids)) {
        if (pos < block.owned_size) {
          update_counts(block.counts, particle_ids.size(), -1);
        }
        erase_in_block(block, pos, size);
        drop_if_empty(it);
        changed();
        return true;
      }
      pos += size;
    }
    return false;
  }

  /**
   * @brief Remove all records of a particle accepted by a predicate.
   * @param pid   Particle that stores the records.
   * @param pred  Callable <tt>bool(BondRecord const &)</tt>.
   * @return number of removed records.
   */
  template <typename Pred> std::size_t remove_if(int pid, Pred &&pred) {
    auto const it = m_blocks.find(pid);
    if (it == m_blocks.end()) {
      return 0u;
    }
    auto &block = it->second;
    std::size_t n_removed = 0u;
    for (std::uint32_t pos = 0u; pos < block.size;) {
      auto const *record = m_data.data() + block.offset + pos;
      auto const size = record_size(record);
      BondRecord const bond(record[0], {record + header_size, size - 2u});
      if (pred(bond)) {
        if (pos < block.owned_size) {
          update_counts(block.counts, bond.particle_ids().size(), -1);
        }
        erase_in_block(block, pos, size);
        ++n_removed;
      } else {
        pos += size;
      }
    }
    if (n_removed) {
      drop_if_empty(it);
      changed();
    }
    return n_removed;
  }

  /** @brief Forget all records of a particle. */
  void erase(int pid) {
    auto const it = m_blocks.find(pid);
    if (it != m_blocks.end()) {
      m_garbage += it->second.capacity;
      m_live -= it->second.size;
      m_blocks.erase(it);
      changed();
    }
  }

  /**
   * @brief Move the records of a particle into a buffer, e.g. to send them
   * along with the particle to another rank. Appends
   * <tt>[size, owned_size, records...]</tt>, also for particles without
   * records, and forgets the records.
   */
  void extract(int pid, std::vector<int> &buffer) {
    auto const *block = find(pid);
    if (block == nullptr) {
      buffer.insert(buffer.end(), {0, 0});
      return;
    }
    buffer.push_back(static_cast<int>(block->size));
    buffer.push_back(static_cast<int>(block->owned_size));
    auto const data = span_of(*block, block->size);
    buffer.insert(buffer.end(), data.begin(), data.end());
    erase(pid);
  }

  /** @brief Append a copy of the records of a particle to a buffer, in the
   *  format of @ref extract(). */
  void copy_to(int pid, std::vector<int> &buffer) const {
    auto const *block = find(pid);
    buffer.push_back(block ? static_cast<int>(block->size) : 0);
    buffer.push_back(block ? static_cast<int>(block->owned_size) : 0);
    if (block) {
      auto const data = span_of(*block, block->size);
      buffer.insert(buffer.end(), data.begin(), data.end());
    }
  }

  /**
   * @brief Store the records of a particle read from a buffer filled by
   * @ref extract(). Replaces the records the particle had.
   * @return the remainder of the buffer.
   */
  std::span<const int> insert(int pid, std::span<const int> buffer) {
    assert(buffer.size() >= header_size);
    auto const size = static_cast<std::uint32_t>(buffer[0]);
    auto const owned_size = static_cast<std::uint32_t>(buffer[1]);
    assert(buffer.size() >= header_size + size and owned_size <= size);
    erase(pid);
    if (size != 0u) {
      auto &block = grow(pid, size);
      std::copy_n(buffer.begin() + header_size, size,
                  m_data.begin() + block.offset);
      block.size = size;
      m_live += size;
      block.owned_size = owned_size;
      block.counts = {};
      for (auto const bond :
           BondRecordRange(span_of(block, block.owned_size))) {
        update_counts(block.counts, bond.particle_ids().size(), +1);
      }
      changed();
    }
    return buffer.subspan(header_size + size);
  }

  /**
   * @brief Rewrite the buffer in a given order of the particles and drop
   * the unused ints.
   *
   * Afterwards, and until the next change, the records of the particle at
   * position @c i of @p pids are available via @ref owned_at() and
   * @ref owned_counts_at(). Particles without records, or not in
   * @p pids, are allowed; the records of the latter are kept.
   *
   * @param pids  Particle ids, in the desired order.
   */
  void compact(std::span<const int> pids) {
    m_ordered.assign(pids.size(), Block{});
    std::vector<int> data;
    data.reserve(m_data.size() - m_garbage);
    for (std::size_t i = 0u; i < pids.size(); ++i) {
      auto const it = m_blocks.find(pids[i]);
      if (it == m_blocks.end()) {
        continue;
      }
      auto &block = it->second;
      auto const offset = static_cast<std::uint32_t>(data.size());
      data.insert(data.end(), m_data.begin() + block.offset,
                  m_data.begin() + block.offset + block.size);
      block.offset = offset;
      block.capacity = block.size;
      m_ordered[i] = block;
      // mark as placed
      block.capacity |= placed_flag;
    }
    // records of particles that are not in the given order
    for (auto &[pid, block] : m_blocks) {
      if (block.capacity & placed_flag) {
        block.capacity &= ~placed_flag;
        continue;
      }
      auto const offset = static_cast<std::uint32_t>(data.size());
      data.insert(data.end(), m_data.begin() + block.offset,
                  m_data.begin() + block.offset + block.size);
      block.offset = offset;
      block.capacity = block.size;
    }
    m_data = std::move(data);
    m_garbage = 0u;
    assert(m_data.size() == m_live);
    m_ordered_valid = true;
  }

  /**
   * @brief Prepare the access by position for @p n particles, without
   * moving any record. Fill it with @ref set_index(), then call
   * @ref index_done().
   */
  void index_resize(std::size_t n) {
    m_ordered.assign(n, Block{});
    m_ordered_valid = false;
  }

  /** @brief Register particle @p pid at position @p i. Safe to call
   *  concurrently for different positions. */
  void set_index(std::size_t i, int pid) {
    auto const *block = find(pid);
    m_ordered[i] = block ? *block : Block{};
  }

  /** @brief Make the access by position available, until the next change. */
  void index_done() { m_ordered_valid = true; }

  /** @brief Whether @ref owned_at() can be used. */
  bool is_indexed() const { return m_ordered_valid; }

  /** @brief Owned records of the particle at position @p i of the last
   *  @ref compact(). */
  BondRecordRange owned_at(std::size_t i) const {
    assert(m_ordered_valid and i < m_ordered.size());
    auto const &block = m_ordered[i];
    return BondRecordRange(span_of(block, block.owned_size));
  }

  /** @brief Owned counts of the particle at position @p i of the last
   *  @ref compact(). */
  Counts const &owned_counts_at(std::size_t i) const {
    assert(m_ordered_valid and i < m_ordered.size());
    return m_ordered[i].counts;
  }

private:
  static constexpr std::uint32_t placed_flag = 1u << 31;

  void drop_if_empty(decltype(m_blocks)::iterator it) {
    if (it->second.size == 0u) {
      m_garbage += it->second.capacity;
      m_blocks.erase(it);
    }
  }
};
