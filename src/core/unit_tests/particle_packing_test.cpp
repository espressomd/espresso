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

#define BOOST_TEST_MODULE "Ghost particle packing test"
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <config/config.hpp>

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleList.hpp"
#include "ghosts.hpp"
#include "ghosts/particle_packing.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <vector>

using namespace GhostComm;

namespace {
/** Two cells of particles, some of which leave the box after the shift. */
struct Fixture {
  BoxGeometry box;
  ParticleList cell1, cell2;
  std::array<ParticleList *, 2> cells{{&cell1, &cell2}};
  Utils::Vector3d const shift{10., 0., -10.};

  Fixture() {
    box.set_length({10., 10., 10.});
    int id = 0;
    for (auto cell : cells) {
      for (int i = 0; i < 3; ++i, ++id) {
        Particle p{};
        p.id() = id;
        p.pos() = {0.5 + 3. * i, 9.5 - 0.7 * id, 0.1 * id};
        p.image_box() = {id, -id, 2 * id};
#ifdef ESPRESSO_BOND_CONSTRAINT
        p.pos_last_time_step() = p.pos() - Utils::Vector3d{0.01, 0.02, 0.03};
#endif
#ifdef ESPRESSO_ROTATION
        p.quat() = Utils::Quaternion<double>{{1., 0.1 * id, 0.2, 0.3}};
#endif
        cell->insert(std::move(p));
      }
    }
  }

  auto pack(unsigned data_parts) {
    CommBuf buf;
    pack_cells(buf, cells, shift, box, data_parts);
    return std::vector<char>(buf.data(), buf.data() + buf.size());
  }

  /** Bytes of the position block of every record, records of @p size. */
  static auto position_blocks(std::vector<char> const &bytes,
                              std::size_t size) {
    auto const block = sizeof(Utils::Vector3d) + sizeof(Utils::Vector3i);
    std::vector<char> out;
    for (std::size_t offset = 0; offset < bytes.size(); offset += size) {
      out.insert(out.end(), bytes.begin() + offset,
                 bytes.begin() + offset + block);
    }
    return out;
  }
};
} // namespace

BOOST_AUTO_TEST_CASE(transmit_size) {
  BoxGeometry box;
  auto const pos_size = sizeof(Utils::Vector3d) + sizeof(Utils::Vector3i);
  BOOST_CHECK_EQUAL(calc_transmit_size(box, GHOSTTRANS_POSITION), pos_size);
  // bonds go to a separate buffer
  BOOST_CHECK_EQUAL(
      calc_transmit_size(box, GHOSTTRANS_POSITION | GHOSTTRANS_BONDS),
      pos_size);
  BOOST_CHECK_EQUAL(
      calc_transmit_size(box, GHOSTTRANS_POSITION | GHOSTTRANS_MOMENTUM),
      calc_transmit_size(box, GHOSTTRANS_POSITION) +
          calc_transmit_size(box, GHOSTTRANS_MOMENTUM));
#ifdef ESPRESSO_BOND_CONSTRAINT
  BOOST_CHECK_EQUAL(
      calc_transmit_size(box, GHOSTTRANS_POSITION | GHOSTTRANS_POS_LAST),
      pos_size + sizeof(Utils::Vector3d));
#endif
  Fixture f;
  BOOST_CHECK_EQUAL(calc_transmit_size(f.cells, f.box, GHOSTTRANS_POSITION),
                    6u * pos_size);
  BOOST_CHECK_EQUAL(calc_transmit_size(f.cells, f.box, GHOSTTRANS_PARTNUM),
                    2u * sizeof(unsigned int));
}

/*
 * pack_cells() packs a plain position push (GHOSTTRANS_POSITION) in a flat
 * loop; with any other data part it uses the generic serializer. Both must
 * write the same position block, since the receiver always decodes with the
 * generic serializer.
 */
BOOST_AUTO_TEST_CASE(flat_position_packer_matches_generic_serializer) {
  Fixture f;
  auto const flat = f.pack(GHOSTTRANS_POSITION);
  BOOST_REQUIRE_EQUAL(flat.size(),
                      6u * calc_transmit_size(f.box, GHOSTTRANS_POSITION));
  // GHOSTTRANS_BONDS adds no bytes to the particle buffer, but takes the
  // generic path
  auto const generic = f.pack(GHOSTTRANS_POSITION | GHOSTTRANS_BONDS);
  BOOST_TEST(flat == generic, boost::test_tools::per_element());
  auto const record = calc_transmit_size(f.box, GHOSTTRANS_POSITION);
  BOOST_TEST(Fixture::position_blocks(flat, record) == flat,
             boost::test_tools::per_element());
#ifdef ESPRESSO_ROTATION
  {
    auto const parts = GHOSTTRANS_POSITION | GHOSTTRANS_QUAT;
    auto const with_quat = f.pack(parts);
    BOOST_TEST(Fixture::position_blocks(
                   with_quat, calc_transmit_size(f.box, parts)) == flat,
               boost::test_tools::per_element());
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  {
    auto const parts = GHOSTTRANS_POSITION | GHOSTTRANS_POS_LAST;
    auto const with_pos_last = f.pack(parts);
    BOOST_TEST(Fixture::position_blocks(
                   with_pos_last, calc_transmit_size(f.box, parts)) == flat,
               boost::test_tools::per_element());
  }
#endif
}

/* The receiver sees the shifted, folded position and image box. */
BOOST_AUTO_TEST_CASE(position_round_trip) {
  Fixture f;
  std::vector<unsigned> const data_parts{
      GHOSTTRANS_POSITION,
#ifdef ESPRESSO_BOND_CONSTRAINT
      GHOSTTRANS_POSITION | GHOSTTRANS_POS_LAST,
#endif
  };
  for (auto const parts : data_parts) {
    CommBuf buf;
    pack_cells(buf, f.cells, f.shift, f.box, parts);
    ParticleList ghost1, ghost2;
    ghost1.resize(f.cell1.size());
    ghost2.resize(f.cell2.size());
    std::array<ParticleList *, 2> ghost_cells{{&ghost1, &ghost2}};
    unpack_cells(buf, ghost_cells, f.box, parts);
    for (std::size_t c = 0; c < 2; ++c) {
      for (std::size_t i = 0; i < f.cells[c]->size(); ++i) {
        auto const &src = f.cells[c]->begin()[i];
        auto const &dst = ghost_cells[c]->begin()[i];
        auto pos = src.pos() + f.shift;
        auto img = src.image_box();
        f.box.fold_position(pos, img);
        BOOST_CHECK_EQUAL(dst.pos(), pos);
        BOOST_CHECK_EQUAL(dst.image_box(), img);
#ifdef ESPRESSO_BOND_CONSTRAINT
        if (parts & GHOSTTRANS_POS_LAST) {
          // sent unshifted
          BOOST_CHECK_EQUAL(dst.pos_last_time_step(), src.pos_last_time_step());
        } else {
          BOOST_CHECK_EQUAL(dst.pos_last_time_step(), Utils::Vector3d{});
        }
#endif
      }
    }
  }
}
