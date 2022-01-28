/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "h5md_core.hpp"

#include "Particle.hpp"
#include "h5md_specification.hpp"
#include "version.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace Writer {
namespace H5md {

using MultiArray3i = boost::multi_array<int, 3>;
using Vector1hs = Utils::Vector<hsize_t, 1>;
using Vector2hs = Utils::Vector<hsize_t, 2>;
using Vector3hs = Utils::Vector<hsize_t, 3>;

static void backup_file(const std::string &from, const std::string &to) {
  /*
   * If the file itself *and* a backup file exists, something must
   * have gone wrong.
   */
  boost::filesystem::path pfrom(from), pto(to);
  try {
    boost::filesystem::copy_file(
        pfrom, pto, boost::filesystem::copy_option::fail_if_exists);
  } catch (const boost::filesystem::filesystem_error &e) {
    throw left_backupfile();
  }
}

template <typename extent_type>
static void extend_dataset(h5xx::dataset &dataset,
                           extent_type const &change_extent) {
  auto const rank = static_cast<h5xx::dataspace>(dataset).rank();
  auto extents = static_cast<h5xx::dataspace>(dataset).extents();
  /* Extend the dataset for another timestep */
  for (int i = 0; i < rank; i++) {
    extents[i] += change_extent[i];
  }
  H5Dset_extent(dataset.hid(), extents.data()); // extend all dims is collective
}

template <typename value_type, typename extent_type>
static void write_dataset(value_type const &data, h5xx::dataset &dataset,
                          extent_type const &change_extent,
                          extent_type const &offset, extent_type const &count) {
  extend_dataset(dataset, change_extent);
  /* write the data to the dataset. */
  h5xx::write_dataset(dataset, data, h5xx::slice(offset, count));
}

void write_script(std::string const &target,
                  boost::filesystem::path const &script_path) {
  std::ifstream scriptfile(script_path.string());
  std::string buffer((std::istreambuf_iterator<char>(scriptfile)),
                     std::istreambuf_iterator<char>());
  auto file = h5xx::file(target, h5xx::file::out);
  auto const group = h5xx::group(file, "parameters/files");
  h5xx::write_attribute(group, "script", buffer);
  file.close();
}

/* Initialize the file related variables after parameters have been set. */
void File::init_file(std::string const &file_path) {
  m_backup_filename = file_path + ".bak";
  boost::filesystem::path script_path(m_script_path);
  m_absolute_script_path = boost::filesystem::canonical(script_path);
  bool file_exists = boost::filesystem::exists(file_path);
  bool backup_file_exists = boost::filesystem::exists(m_backup_filename);
  /* Perform a barrier synchronization. Otherwise one process might already
   * create the file while another still checks for its existence. */
  m_comm.barrier();
  if (file_exists) {
    if (H5MD_Specification::is_compliant(file_path)) {
      /*
       * If the file exists and has a valid H5MD structure, let's create a
       * backup of it. This has the advantage, that the new file can
       * just be deleted if the simulation crashes at some point and we
       * still have a valid trajectory, we can start from.
       */
      if (m_comm.rank() == 0)
        backup_file(file_path, m_backup_filename);
      load_file(file_path);
    } else {
      throw incompatible_h5mdfile();
    }
  } else {
    if (backup_file_exists)
      throw left_backupfile();
    create_file(file_path);
  }
}

void File::load_datasets() {
  for (auto const &d : H5MD_Specification::DATASETS) {
    if (d.is_link)
      continue;
    datasets[d.path()] = h5xx::dataset(m_h5md_file, d.path());
  }
}

void File::create_groups() {
  h5xx::group group(m_h5md_file);
  for (auto const &d : H5MD_Specification::DATASETS) {
    h5xx::group new_group(group, d.group);
  }
}

static std::vector<hsize_t> create_dims(hsize_t rank, hsize_t data_dim) {
  switch (rank) {
  case 3:
    return std::vector<hsize_t>{0, 0, data_dim};
  case 2:
    return std::vector<hsize_t>{0, data_dim};
  case 1:
    return std::vector<hsize_t>{data_dim};
  default:
    throw std::runtime_error(
        "H5MD Error: datasets with this dimension are not implemented\n");
  }
}

static std::vector<hsize_t> create_chunk_dims(hsize_t rank, hsize_t data_dim) {
  hsize_t chunk_size = (rank > 1) ? 1000 : 1;
  switch (rank) {
  case 3:
    return {1, chunk_size, data_dim};
  case 2:
    return {1, chunk_size};
  case 1:
    return {chunk_size};
  default:
    throw std::runtime_error(
        "H5MD Error: datasets with this dimension are not implemented\n");
  }
}

void File::create_datasets() {
  namespace hps = h5xx::policy::storage;
  for (const auto &d : H5MD_Specification::DATASETS) {
    if (d.is_link)
      continue;
    auto maxdims = std::vector<hsize_t>(d.rank, H5S_UNLIMITED);
    auto dataspace = h5xx::dataspace(create_dims(d.rank, d.data_dim), maxdims);
    auto storage = hps::chunked(create_chunk_dims(d.rank, d.data_dim))
                       .set(hps::fill_value(-10));
    datasets[d.path()] = h5xx::dataset(m_h5md_file, d.path(), d.type, dataspace,
                                       storage, H5P_DEFAULT, H5P_DEFAULT);
  }
}

void File::load_file(const std::string &file_path) {
  m_h5md_file = h5xx::file(file_path, m_comm, MPI_INFO_NULL, h5xx::file::out);
  load_datasets();
}

static void write_box(const BoxGeometry &geometry, h5xx::dataset &dataset) {
  auto const extents = static_cast<h5xx::dataspace>(dataset).extents();
  extend_dataset(dataset, Vector2hs{1, 0});
  h5xx::write_dataset(dataset, geometry.length(),
                      h5xx::slice(Vector2hs{extents[0], 0}, Vector2hs{1, 3}));
}

void write_attributes(const std::string &espresso_version,
                      h5xx::file &h5md_file) {
  auto h5md_group = h5xx::group(h5md_file, "h5md");
  h5xx::write_attribute(h5md_group, "version",
                        boost::array<hsize_t, 2>{{1, 1}});
  auto h5md_creator_group = h5xx::group(h5md_group, "creator");
  h5xx::write_attribute(h5md_creator_group, "name", "ESPResSo");
  h5xx::write_attribute(h5md_creator_group, "version", espresso_version);
  auto h5md_author_group = h5xx::group(h5md_group, "author");
  h5xx::write_attribute(h5md_author_group, "name", "N/A");
  auto group = h5xx::group(h5md_file, "particles/atoms/box");
  h5xx::write_attribute(group, "dimension", 3);
  h5xx::write_attribute(group, "boundary", "periodic");
}

void File::write_units() {
  h5xx::write_attribute(datasets["particles/atoms/mass/value"], "unit",
                        m_mass_unit);
  h5xx::write_attribute(datasets["particles/atoms/charge/value"], "unit",
                        m_charge_unit);
  h5xx::write_attribute(datasets["particles/atoms/position/value"], "unit",
                        m_length_unit);
  h5xx::write_attribute(datasets["particles/atoms/velocity/value"], "unit",
                        m_velocity_unit);
  h5xx::write_attribute(datasets["particles/atoms/force/value"], "unit",
                        m_force_unit);
  h5xx::write_attribute(datasets["particles/atoms/id/time"], "unit",
                        m_time_unit);
}

void hard_link(h5xx::file const &file, std::string from, std::string to) {
  if (H5Lcreate_hard(file.hid(), from.c_str(), file.hid(), to.c_str(),
                     H5P_DEFAULT, H5P_DEFAULT) < 0) {
    throw std::runtime_error("Error creating hard link for " + to);
  }
}

void File::create_hard_links() {
  std::string path_step = "particles/atoms/id/step";
  std::string path_time = "particles/atoms/id/time";
  for (auto &ds : H5MD_Specification::DATASETS) {
    if (ds.name == "step" and ds.is_link) {
      hard_link(m_h5md_file, path_step, ds.path());
    } else if (ds.name == "time" and ds.is_link) {
      hard_link(m_h5md_file, path_time, ds.path());
    }
  }
}

void File::create_file(const std::string &file_path) {
  if (m_comm.rank() == 0)
    write_script(file_path, m_absolute_script_path);
  m_comm.barrier();
  m_h5md_file = h5xx::file(file_path, m_comm, MPI_INFO_NULL, h5xx::file::out);
  create_groups();
  create_datasets();
  write_attributes(ESPRESSO_VERSION, m_h5md_file);
  write_units();
  create_hard_links();
}

void File::close() {
  if (m_comm.rank() == 0)
    boost::filesystem::remove(m_backup_filename);
}

namespace detail {

template <std::size_t rank> struct slice_info {};

template <> struct slice_info<3> {
  static auto extent(hsize_t n_part_diff) {
    return Vector3hs{1, n_part_diff, 0};
  }
  static constexpr auto count() { return Vector3hs{1, 1, 3}; }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector3hs{n_time_steps, prefix, 0};
  }
};

template <> struct slice_info<2> {
  static auto extent(hsize_t n_part_diff) { return Vector2hs{1, n_part_diff}; }
  static constexpr auto count() { return Vector2hs{1, 1}; }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector2hs{n_time_steps, prefix};
  }
};

} // namespace detail
template <std::size_t dim, typename Op>
void write_td_particle_property(hsize_t prefix, hsize_t n_part_global,
                                ParticleRange const &particles,
                                h5xx::dataset &dataset, Op op) {
  auto const old_extents = static_cast<h5xx::dataspace>(dataset).extents();
  auto const extent_particle_number =
      std::max(n_part_global, old_extents[1]) - old_extents[1];
  extend_dataset(dataset,
                 detail::slice_info<dim>::extent(extent_particle_number));
  auto const count = detail::slice_info<dim>::count();
  auto offset = detail::slice_info<dim>::offset(old_extents[0], prefix);
  for (auto const &p : particles) {
    h5xx::write_dataset(dataset, op(p), h5xx::slice(offset, count));
    // advance in the particle dimension
    offset[1] += 1;
  }
}

void File::write(const ParticleRange &particles, double time, int step,
                 BoxGeometry const &geometry) {
  write_box(geometry, datasets["particles/atoms/box/edges/value"]);
  write_connectivity(particles);

  int const n_part_local = particles.size();
  // calculate count and offset
  int prefix = 0;
  // calculate prefix for write of the current process
  BOOST_MPI_CHECK_RESULT(MPI_Exscan,
                         (&n_part_local, &prefix, 1, MPI_INT, MPI_SUM, m_comm));
  auto const extents =
      static_cast<h5xx::dataspace>(datasets["particles/atoms/id/value"])
          .extents();

  auto const n_part_global =
      boost::mpi::all_reduce(m_comm, n_part_local, std::plus<int>());

  write_td_particle_property<2>(
      prefix, n_part_global, particles, datasets["particles/atoms/id/value"],
      [](auto const &p) { return Utils::Vector<int, 1>{p.p.identity}; });
  write_dataset(Utils::Vector<double, 1>{time},
                datasets["particles/atoms/id/time"], Vector1hs{1},

                Vector1hs{extents[0]}, Vector1hs{1});
  write_dataset(Utils::Vector<int, 1>{step},
                datasets["particles/atoms/id/step"], Vector1hs{1},
                Vector1hs{extents[0]}, Vector1hs{1});

  write_td_particle_property<2>(
      prefix, n_part_global, particles,
      datasets["particles/atoms/species/value"],
      [](auto const &p) { return Utils::Vector<int, 1>{p.p.type}; });

  write_td_particle_property<2>(
      prefix, n_part_global, particles, datasets["particles/atoms/mass/value"],
      [](auto const &p) { return Utils::Vector<double, 1>{p.p.mass}; });

  write_td_particle_property<3>(
      prefix, n_part_global, particles,
      datasets["particles/atoms/position/value"],
      [&](auto const &p) { return folded_position(p.r.p, geometry); });
  write_td_particle_property<3>(prefix, n_part_global, particles,
                                datasets["particles/atoms/image/value"],
                                [](auto const &p) { return p.l.i; });

  write_td_particle_property<3>(prefix, n_part_global, particles,
                                datasets["particles/atoms/velocity/value"],
                                [](auto const &p) { return p.m.v; });

  write_td_particle_property<3>(prefix, n_part_global, particles,
                                datasets["particles/atoms/force/value"],
                                [](auto const &p) { return p.f.f; });
  write_td_particle_property<2>(
      prefix, n_part_global, particles,
      datasets["particles/atoms/charge/value"],
      [](auto const &p) { return Utils::Vector<double, 1>{p.p.q}; });
}
void File::write_connectivity(const ParticleRange &particles) {
  MultiArray3i bond(boost::extents[0][0][0]);
  int particle_index = 0;
  for (auto const &p : particles) {
    int nbonds_local = bond.shape()[1];
    for (auto const b : p.bonds()) {
      auto const partner_ids = b.partner_ids();
      if (partner_ids.size() == 1) {
        bond.resize(boost::extents[1][nbonds_local + 1][2]);
        bond[0][nbonds_local][0] = p.p.identity;
        bond[0][nbonds_local][1] = partner_ids[0];
        nbonds_local++;
      }
    }
    particle_index++;
  }

  int n_bonds_local = bond.shape()[1];
  int prefix_bonds = 0;
  BOOST_MPI_CHECK_RESULT(
      MPI_Exscan, (&n_bonds_local, &prefix_bonds, 1, MPI_INT, MPI_SUM, m_comm));
  auto const n_bonds_total =
      boost::mpi::all_reduce(m_comm, n_bonds_local, std::plus<int>());
  auto const extents =
      static_cast<h5xx::dataspace>(datasets["connectivity/atoms/value"])
          .extents();
  Vector3hs offset_bonds = {extents[0], static_cast<hsize_t>(prefix_bonds), 0};
  Vector3hs count_bonds = {1, static_cast<hsize_t>(n_bonds_local), 2};
  auto const n_bond_diff =
      std::max(static_cast<hsize_t>(n_bonds_total), extents[1]) - extents[1];
  Vector3hs change_extent_bonds = {1, static_cast<hsize_t>(n_bond_diff), 0};
  write_dataset(bond, datasets["connectivity/atoms/value"], change_extent_bonds,
                offset_bonds, count_bonds);
}

void File::flush() { m_h5md_file.flush(); }

} /* namespace H5md */
} /* namespace Writer */
