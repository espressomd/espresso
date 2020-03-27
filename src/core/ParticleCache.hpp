/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef CORE_PARTICLE_CACHE_HPP
#define CORE_PARTICLE_CACHE_HPP

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include <boost/version.hpp>
/* Work around bug in boost, see
   https://github.com/boostorg/container/commit/5e4a107e82ab3281688311d22d2bfc2fddcf84a3
*/
#if BOOST_VERSION < 106400
#include <boost/container/detail/pair.hpp>
#endif

#include <boost/container/flat_set.hpp>
#include <boost/mpi/collectives.hpp>

#include "MpiCallbacks.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include <boost/algorithm/clamp.hpp>
#include <utils/NoOp.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/flat_set.hpp>

namespace detail {
/**
 * @brief Compare particles by id.
 */
class IdCompare {
public:
  template <typename Particle>
  bool operator()(Particle const &a, Particle const &b) const {
    return a.identity() < b.identity();
  }
};

/**
 * @brief Merge two ordered containers into a new one.
 *
 * This implementation has a different tradeoff than
 * flat_set::merge, here we use O(N) extra memory to
 * get O(N) time complexity, while the flat_map implementation
 * avoids extra memory usage, but will cause O(N^2) copies
 * on average.
 * Inspired by the example implementation in
 * http://en.cppreference.com/w/cpp/algorithm/merge.
 */
template <typename Container, typename Compare> class Merge {
  Compare m_comp;

public:
  explicit Merge(Compare &&comp = Compare{}) : m_comp(comp) {}
  Container operator()(Container const &a, Container const &b) const {
    Container ret;
    ret.reserve(a.size() + b.size());

    auto first1 = a.begin();
    auto last1 = a.end();
    auto first2 = b.begin();
    auto last2 = b.end();

    while (first1 != last1) {
      if (first2 == last2) {
        /* The 2nd range has no more elements, so we can
           just copy the rest of range 1. */
        for (; first1 != last1; ++first1) {
          ret.emplace_hint(ret.end(), *first1);
        }
        break;
      }

      if (m_comp(*first2, *first1)) {
        ret.emplace_hint(ret.end(), *first2);
        ++first2;
      } else {
        ret.emplace_hint(ret.end(), *first1);
        ++first1;
      }
    }

    /* The first range has no more elements, so we can
       just copy the rest of range 2. */
    for (; first2 != last2; ++first2)
      ret.emplace_hint(ret.end(), *first2);

    return ret;
  }
};
} // namespace detail

/* Mark merge as commutative with all containers */
namespace boost {
namespace mpi {
template <typename Container>
struct is_commutative<::detail::Merge<Container, ::detail::IdCompare>,
                      Container> : public boost::mpl::true_ {};
} // namespace mpi
} // namespace boost

/**
 * @brief Particle cache on the master.
 *
 * This class implements cached access to all particles in a
 * particle range on the master node.
 * This implementation fetches all particles to
 * the master on first access. Updates of the particle data are
 * triggered automatically on access. The data in the cache
 * is invalidated automatically on_particle_change, and then
 * updated on the next access.
 *
 * To update the cache particles are sorted by id on the nodes,
 * and the sorted arrays a merged in a reduction tree, until the
 * master node receives a complete and sorted particle array.
 *
 * This class can be customized by running a unary operation on
 * the particles. This op is run on all the nodes. It can be used
 * e.g. to fold or unfold the coordinates on the fly.
 *
 * To iterate over the particles using the iterators is more
 * efficient than using operator[].
 *
 * All functions in the public interface can only be called on
 * the master node.
 */
template <typename GetParticles, typename UnaryOp = Utils::NoOp,
          typename Range = typename std::remove_reference<decltype(
              std::declval<GetParticles>()())>::type,
          typename Particle = typename std::iterator_traits<
              typename Range::iterator>::value_type>
class ParticleCache {
  using map_type = std::vector<Particle>;
  /* Callback system we're on */
  Communication::MpiCallbacks &m_cb;

  /** Index mapping particle ids to the position
      in remote_parts. */
  std::unordered_map<int, int> id_index;
  /** The particle data */
  std::vector<Particle> remote_parts;
  /** State */
  bool m_valid;

  Communication::CallbackHandle<> update_cb;

  /** Functor to get a particle range */
  GetParticles m_parts;
  /** Functor which is applied to the
      particles before they are gathered,
      e.g. position folding */
  UnaryOp m_op;

  /**
   * @brief Actual update implementation.
   *
   * This gets a new particle range, packs
   * the particles into a buffer and then
   * merges these buffers hierarchically to the
   * master node
   */
  void m_update() {}

  void m_update_index() {
    /* Try to avoid rehashing along the way */
    id_index.reserve(remote_parts.size() + 1);

    int index = 0;
    for (auto const &p : remote_parts) {
      id_index.insert(std::make_pair(p.identity(), index++));
    }
  }

public:
  using value_iterator = typename map_type::const_iterator;
  using value_type = Particle;

  ParticleCache() = delete;
  ParticleCache(Communication::MpiCallbacks &cb, GetParticles parts,
                UnaryOp &&op = UnaryOp{})
      : m_cb(cb), m_valid(false), update_cb(&cb, [this]() { m_update(); }),
        m_parts(parts), m_op(std::forward<UnaryOp>(op)) {}
  /* Because the this ptr is captured by the callback lambdas,
   * this class can be neither copied nor moved. */
  ParticleCache(ParticleCache const &) = delete;
  ParticleCache operator=(ParticleCache const &) = delete;
  ParticleCache(ParticleCache &&) = delete;
  ParticleCache operator=(ParticleCache &&) = delete;

  /**
   * @brief Clear cache.
   */
  void clear() {
    id_index.clear();
    remote_parts.clear();
  }

  /**
   * @brief Iterator pointing to the particle with the lowest
   * id.
   *
   * Returns a random access iterator that traverses the
   * particles
   * in order of ascending id. If the cache is not up-to-date,
   * an update is triggered. This iterator stays valid as long
   * as the cache is valid. Since the cache could be invalidated
   * and updated elsewhere, iterators into the cache should not
   * be stored.
   */
  value_iterator begin() {
    assert(m_cb.comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.begin());
  }

  /**
   * @brief Iterator pointing past the particle with the highest
   * id.
   *
   * If the cache is not up-to-date,
   * an update is triggered.
   */
  value_iterator end() {
    assert(m_cb.comm().rank() == 0);

    if (!m_valid)
      update();

    return value_iterator(remote_parts.end());
  }

  /**
   * @brief Returns true if the cache is up-to-date.
   *
   * If false, particle access will trigger an update.
   */
  bool valid() const { return m_valid; }

  /**
   * @brief Invalidate the cache and free memory.
   */
  void invalidate() {
    clear();
    /* Release memory */
    remote_parts.shrink_to_fit();
    /* Adjust state */
    m_valid = false;
  }

  /**
   * @brief Update particle information.
   *
   * This triggers a global update. All nodes
   * sort their particle by id, and send them
   * to the master.
   *
   * Complexity: 2*N*(1 - 0.5^(log(p) + 1))
   */
  void update() {
    if (m_valid)
      return;

    remote_parts.clear();

    auto const ids = get_particle_ids();
    auto const chunk_size = fetch_cache_max_size();

    for (size_t offset = 0; offset < ids.size();) {
      auto const this_size =
          boost::algorithm::clamp(chunk_size, 0, ids.size() - offset);
      auto const chunk_ids =
          Utils::make_const_span(ids.data() + offset, this_size);

      prefetch_particle_data(chunk_ids);

      for (auto id : chunk_ids) {
        remote_parts.push_back(get_particle_data(id));

        auto &p = remote_parts.back();
        p.r.p += image_shift(p.l.i, box_geo.length());
        p.l.i = {};
      }

      offset += this_size;
    }

    m_update_index();

    m_valid = true;
  }

  /** Number of particles in the config.
   *
   * Complexity: O(1)
   */
  size_t size() {
    assert(m_cb.comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.size();
  }

  /**
   * @brief size() == 0 ?
   *
   * Complexity: O(1)
   */
  bool empty() {
    assert(m_cb.comm().rank() == 0);

    if (!m_valid)
      update();

    return remote_parts.empty();
  }

  /**
   * @brief Access particle by id.
   * If the particle config is not valid this will trigger
   * a global update.
   * Will throw std::out_of_range if the particle does
   * not exists.
   *
   * Complexity: O(1)
   */
  Particle const &operator[](int id) {
    assert(m_cb.comm().rank() == 0);

    if (!m_valid)
      update();

    /* Fetch the particle using the index. Here
       begin()[] with the position is used to
       get constant time access. remote_parts[id]
       would also be correct, but is O(n*log(n)). */
    return remote_parts.begin()[id_index.at(id)];
  }
};

#endif
