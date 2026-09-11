/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
/** \file
 *  Reusable particle packing/unpacking for ghost communications.
 *
 *  A cell's particles are @ref Particle views over a contiguous range of
 *  @ref ParticleStore rows, so the loops below iterate @ref Cell::particles
 *  and read/write the store columns through the view accessors. Ghost-ness is
 *  structural (a row's position in the store), so the PARTNUM bootstrap resizes
 *  a ghost cell through @ref CellParticleStorage::resize_ghost_storage instead
 *  of flagging individual particles.
 */

#include "particle_packing.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/Cell.hpp"
#include "cell_system/ParticleListOperations.hpp"
#include "ghosts.hpp"
#include "particle_store/ParticleStore.hpp"

#include <utils/Vector.hpp>
#include <utils/serialization/memcpy_archive.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

namespace GhostComm {

/** @brief Type of reduction to carry out during serialization. */
enum class ReductionPolicy {
  /** @brief Reduction for domain-to-domain particle communication. */
  MOVE,
  /** @brief Reduction for cell-to-cell particle update. */
  UPDATE,
};

/** @brief Whether to save the state to or load the state from the archive. */
enum class SerializationDirection { SAVE, LOAD };

/**
 * @brief Serialize particle data, possibly with reduction.
 * The reduction can take place during the save stage, e.g. to apply
 * a ghost shift to the particle position, or during the load stage,
 * e.g. to transfer momentum between particles in two local cells.
 *
 * Every field goes through the @ref Particle view accessors, which read and
 * write the @ref ParticleStore columns of the view's row. The scalar accessors
 * return references and can be archived directly; the vector accessors return
 * proxies, so those fields are materialized into a local value first and
 * written back explicitly on load. The wire layout is identical either way.
 */
template <class Archive>
static void
serialize_and_reduce(Archive &ar, Particle &p, unsigned int data_parts,
                     ReductionPolicy policy, SerializationDirection direction,
                     BoxGeometry const &box_geo,
                     Utils::Vector3d const *ghost_shift) {
  auto const loading = direction == SerializationDirection::LOAD;

  /** Archive a field that the accessor exposes as a proxy rather than a
   *  reference: materialize on save, write back on load. The value type comes
   *  from the CONST accessor overload (which returns by value); the non-const
   *  overload hands back the writable column proxy. */
  auto archive_vector = [&ar, &p, loading](auto &&accessor) {
    using value_type =
        std::remove_cvref_t<decltype(accessor(std::as_const(p)))>;
    value_type value{};
    if (not loading) {
      value = accessor(p);
    }
    ar & value;
    if (loading) {
      accessor(p) = value;
    }
  };

  if (data_parts & GHOSTTRANS_PROPRTS) {
    ar & p.id() & p.mol_id() & p.type() & p.propagation();
#ifdef ESPRESSO_ROTATION
    ar & p.rotation();
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    archive_vector([](auto &&q) -> decltype(auto) { return q.rinertia(); });
#endif
#endif
#ifdef ESPRESSO_MASS
    ar & p.mass();
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    ar & p.q();
#endif
#ifdef ESPRESSO_DIPOLES
    ar & p.dipm();
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    archive_vector([](auto &&q) -> decltype(auto) { return q.mu_E(); });
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    ar & p.vs_relative();
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    archive_vector([](auto &&q) -> decltype(auto) { return q.gamma(); });
#ifdef ESPRESSO_ROTATION
    archive_vector([](auto &&q) -> decltype(auto) { return q.gamma_rot(); });
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    ar & p.fixed();
    archive_vector([](auto &&q) -> decltype(auto) { return q.ext_force(); });
#ifdef ESPRESSO_ROTATION
    archive_vector([](auto &&q) -> decltype(auto) { return q.ext_torque(); });
#endif
#endif
#ifdef ESPRESSO_ENGINE
    ar & p.swimming();
#endif
  }
  if (data_parts & GHOSTTRANS_POSITION) {
    if (not loading and ghost_shift != nullptr) {
      auto pos = Utils::Vector3d(p.pos()) + *ghost_shift;
      auto img = Utils::Vector3i(p.image_box());
      box_geo.fold_position(pos, img);
      ar & pos;
      ar & img;
    } else {
      archive_vector([](auto &&q) -> decltype(auto) { return q.pos(); });
      archive_vector([](auto &&q) -> decltype(auto) { return q.image_box(); });
    }
#ifdef ESPRESSO_BOND_CONSTRAINT
    archive_vector(
        [](auto &&q) -> decltype(auto) { return q.pos_last_time_step(); });
#endif
  }
  // Wire-symmetry: pack and unpack use the same data_parts per exchange,
  // so the layout is always symmetric; QUAT follows POSITION in the stream
  // when both are set (position push), and TORQUE follows FORCE when both
  // are set (force reduce).
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_QUAT) {
    archive_vector([](auto &&q) -> decltype(auto) { return q.quat(); });
  }
#endif
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    archive_vector([](auto &&q) -> decltype(auto) { return q.v(); });
#ifdef ESPRESSO_ROTATION
    archive_vector([](auto &&q) -> decltype(auto) { return q.omega(); });
#endif
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    if (policy == ReductionPolicy::UPDATE and loading) {
      Utils::Vector3d force;
      ar & force;
      p.force() += force;
    } else {
      archive_vector([](auto &&q) -> decltype(auto) { return q.force(); });
    }
  }
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_TORQUE) {
    if (policy == ReductionPolicy::UPDATE and loading) {
      Utils::Vector3d torque;
      ar & torque;
      p.torque() += torque;
    } else {
      archive_vector([](auto &&q) -> decltype(auto) { return q.torque(); });
    }
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (data_parts & GHOSTTRANS_RATTLE) {
    if (policy == ReductionPolicy::UPDATE and loading) {
      Utils::Vector3d correction;
      ar & correction;
      p.rattle_correction() += correction;
    } else {
      archive_vector(
          [](auto &&q) -> decltype(auto) { return q.rattle_correction(); });
    }
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (data_parts & GHOSTTRANS_DIPFLD) {
    if (policy == ReductionPolicy::UPDATE and loading) {
      Utils::Vector3d dip_fld;
      ar & dip_fld;
      p.dip_fld() += dip_fld;
    } else {
      archive_vector([](auto &&q) -> decltype(auto) { return q.dip_fld(); });
    }
  }
#endif
}

std::size_t calc_transmit_size(BoxGeometry const &, unsigned data_parts) {
  /* Compositional, compile-time-constant size of the wire layout, matching the
   * field order and ifdef structure of serialize_and_reduce exactly. Each
   * `ar & <field>` packs sizeof(field) tightly (the MemcpyArchive uses memcpy
   * for trivially/bitwise-serializable types, with no alignment padding between
   * fields), so summing sizeof(T) reproduces the packed size byte for byte. The
   * cold PODs are bitwise-serializable and are packed whole, INCLUDING internal
   * struct padding (e.g. VirtualSitesRelativeParameters is 80 B: int + 4 B pad
   * + double + 2 quaternions, NOT 76).
   *
   * A live @ref Particle cannot serve as a sizing dummy here the way it does
   * for the AoS layout: a default-constructed view is not attached to a store,
   * so every accessor would assert. */
  std::size_t size = 0ul;
  if (data_parts & GHOSTTRANS_PROPRTS) {
    size += 4ul * sizeof(int); // id, mol_id, type, propagation
#ifdef ESPRESSO_ROTATION
    size += sizeof(std::uint8_t); // rotation
#ifdef ESPRESSO_ROTATIONAL_INERTIA
    size += sizeof(Utils::Vector3d); // rinertia
#endif
#endif
#ifdef ESPRESSO_MASS
    size += sizeof(double); // mass
#endif
#ifdef ESPRESSO_ELECTROSTATICS
    size += sizeof(double); // q
#endif
#ifdef ESPRESSO_DIPOLES
    size += sizeof(double); // dipm
#endif
#ifdef ESPRESSO_LB_ELECTROHYDRODYNAMICS
    size += sizeof(Utils::Vector3d); // mu_E
#endif
#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
    size += sizeof(VirtualSitesRelativeParameters); // vs_relative
#endif
#ifdef ESPRESSO_THERMOSTAT_PER_PARTICLE
    size += sizeof(ParticleStore::GammaValue); // gamma
#ifdef ESPRESSO_ROTATION
    size += sizeof(ParticleStore::GammaValue); // gamma_rot
#endif
#endif
#ifdef ESPRESSO_EXTERNAL_FORCES
    size += sizeof(std::uint8_t);    // ext_flag / fixed
    size += sizeof(Utils::Vector3d); // ext_force
#ifdef ESPRESSO_ROTATION
    size += sizeof(Utils::Vector3d); // ext_torque
#endif
#endif
#ifdef ESPRESSO_ENGINE
    size += sizeof(ParticleParametersSwimming); // swimming
#endif
  }
  if (data_parts & GHOSTTRANS_POSITION) {
    size += sizeof(Utils::Vector3d); // pos
    size += sizeof(Utils::Vector3i); // image_box
#ifdef ESPRESSO_BOND_CONSTRAINT
    size += sizeof(Utils::Vector3d); // pos_last_time_step
#endif
  }
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_QUAT) {
    size += sizeof(Utils::Quaternion<double>); // quat
  }
#endif
  if (data_parts & GHOSTTRANS_MOMENTUM) {
    size += sizeof(Utils::Vector3d); // v
#ifdef ESPRESSO_ROTATION
    size += sizeof(Utils::Vector3d); // omega
#endif
  }
  if (data_parts & GHOSTTRANS_FORCE) {
    size += sizeof(Utils::Vector3d); // force
  }
#ifdef ESPRESSO_ROTATION
  if (data_parts & GHOSTTRANS_TORQUE) {
    size += sizeof(Utils::Vector3d); // torque
  }
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  if (data_parts & GHOSTTRANS_RATTLE) {
    size += sizeof(Utils::Vector3d); // rattle_correction
  }
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  if (data_parts & GHOSTTRANS_DIPFLD) {
    size += sizeof(Utils::Vector3d); // dip_fld
  }
#endif
  /* GHOSTTRANS_BONDS contributes nothing here: the bond payload travels in the
   * dedicated ragged bond buffer (CommBuf::bondbuf), sized at pack time. */
  return size;
}

std::size_t calc_transmit_size(std::span<Cell *const> cells,
                               BoxGeometry const &box_geo,
                               unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM)
    return sizeof(unsigned int) * cells.size();

  // Count the rows the pack/unpack loops actually visit: @ref Cell::particles
  // hands out views over the COMMITTED rows only (skipping pending-removed
  // ones), whereas @ref Cell::size also counts staged, not-yet-committed
  // particles -- which is what the PARTNUM bootstrap above must report, but
  // would oversize the data buffer here.
  auto const n_part =
      std::accumulate(cells.begin(), cells.end(), std::size_t{0},
                      [](std::size_t sum, auto const *cell) {
                        return sum + cell->rows().size();
                      });

  return n_part * calc_transmit_size(box_geo, data_parts);
}

void pack_cells(CommBuf &buf, std::span<Cell *const> cells,
                Utils::Vector3d const &shift, BoxGeometry const &box_geo,
                unsigned data_parts) {
  /* reallocate send buffer */
  buf.resize(calc_transmit_size(cells, box_geo, data_parts));
  buf.bonds().clear();

  auto archiver = Utils::MemcpyOArchive{buf.make_span()};

  /* Construct archive that pushes back to the bond buffer */
  namespace io = boost::iostreams;
  io::stream<io::back_insert_device<std::vector<char>>> os{
      io::back_inserter(buf.bonds())};
  boost::archive::binary_oarchive bond_archiver{os};

  /* put in data */
  for (auto *cell : cells) {
    if (data_parts & GHOSTTRANS_PARTNUM) {
      // Total count (committed rows + staged): a source cell resized earlier in
      // this pass as a ghost destination has its ghosts staged, not yet
      // committed. The receiver resizes to this pending count.
      assert(cell->size() <= std::numeric_limits<unsigned int>::max());
      auto np = static_cast<unsigned int>(cell->size());
      archiver << np;
    } else {
      for (auto &p : cell->particles()) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::SAVE, box_geo, &shift);
        if (data_parts & GHOSTTRANS_BONDS) {
          bond_archiver << p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_written() == buf.size());
}

void unpack_cells(CommBuf &buf, std::span<Cell *const> cells,
                  BoxGeometry const &box_geo, unsigned data_parts) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};

  if (data_parts & GHOSTTRANS_PARTNUM) {
    for (auto *cell : cells) {
      unsigned int np;
      archiver >> np;
      // Ghost-ness is structural (the row's position in the store), so the
      // resize both sizes the cell and makes its rows ghost rows.
      CellParticleStorage::resize_ghost_storage(*cell, np);
    }
  } else {
    for (auto *cell : cells) {
      for (auto &p : cell->particles()) {
        serialize_and_reduce(archiver, p, data_parts, ReductionPolicy::MOVE,
                             SerializationDirection::LOAD, box_geo, nullptr);
      }
    }
    if (data_parts & GHOSTTRANS_BONDS) {
      namespace io = boost::iostreams;
      io::stream<io::array_source> bond_stream(
          io::array_source{buf.bonds().data(), buf.bonds().size()});
      boost::archive::binary_iarchive bond_archiver(bond_stream);

      for (auto *cell : cells) {
        for (auto &p : cell->particles()) {
          bond_archiver >> p.bonds();
        }
      }
    }
  }

  assert(archiver.bytes_read() == buf.size());

  buf.bonds().clear();
}

void add_forces(CommBuf &buf, std::span<Cell *const> cells,
                unsigned data_parts) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto *cell : cells) {
    for (auto &part : cell->particles()) {
      if (data_parts & GHOSTTRANS_FORCE) {
        Utils::Vector3d force;
        archiver >> force;
        part.force() += force;
      }
#ifdef ESPRESSO_ROTATION
      if (data_parts & GHOSTTRANS_TORQUE) {
        Utils::Vector3d torque;
        archiver >> torque;
        part.torque() += torque;
      }
#endif
    }
  }
}

#ifdef ESPRESSO_BOND_CONSTRAINT
void add_rattle(CommBuf &buf, std::span<Cell *const> cells) {
  /* put back data */
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto *cell : cells) {
    for (auto &part : cell->particles()) {
      Utils::Vector3d correction;
      archiver >> correction;
      part.rattle_correction() += correction;
    }
  }
}
#endif

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
void add_dip_fld(CommBuf &buf, std::span<Cell *const> cells) {
  auto archiver = Utils::MemcpyIArchive{buf.make_span()};
  for (auto *cell : cells) {
    for (auto &part : cell->particles()) {
      Utils::Vector3d dip_fld;
      archiver >> dip_fld;
      part.dip_fld() += dip_fld;
    }
  }
}
#endif

void local_cell_copy(Cell &src, Cell &dst, Utils::Vector3d const &shift,
                     BoxGeometry const &box_geo, unsigned data_parts) {
  if (data_parts & GHOSTTRANS_PARTNUM) {
    CellParticleStorage::resize_ghost_storage(dst, src.size());
    return;
  }

  // Committed rows only: the lockstep loop below walks Cell::particles.
  assert(src.rows().size() == dst.rows().size());

  CommBuf buffer;
  buffer.resize(calc_transmit_size(box_geo, data_parts));

  auto src_range = src.particles();
  auto dst_range = dst.particles();
  auto src_it = src_range.begin();
  auto dst_it = dst_range.begin();
  for (; src_it != src_range.end() and dst_it != dst_range.end();
       ++src_it, ++dst_it) {
    auto ar_out = Utils::MemcpyOArchive{buffer.make_span()};
    auto ar_in = Utils::MemcpyIArchive{buffer.make_span()};
    auto &p1 = *src_it;
    auto &p2 = *dst_it;
    serialize_and_reduce(ar_out, p1, data_parts, ReductionPolicy::UPDATE,
                         SerializationDirection::SAVE, box_geo, &shift);
    serialize_and_reduce(ar_in, p2, data_parts, ReductionPolicy::UPDATE,
                         SerializationDirection::LOAD, box_geo, nullptr);
    if (data_parts & GHOSTTRANS_BONDS) {
      p2.bonds() = p1.bonds();
    }
  }
}

} // namespace GhostComm
