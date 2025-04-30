/*
 i Copyright (C) 2010-2022 The ESPResSo project
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

#include "hdf5_patches.hpp" // must appear first

#include "h5md_core.hpp"
#include "h5md_dataset.hpp"
#include "h5md_specification.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "lees_edwards/LeesEdwardsBC.hpp"

#include "config/version.hpp"

#include <utils/Vector.hpp>

#include <boost/array.hpp>
#include <boost/filesystem.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/multi_array.hpp>

#include <h5xx/h5xx.hpp>
#include <highfive/highfive.hpp>
#define H5_USE_BOOST
#include <highfive/boost.hpp>

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace h5xx {
template <typename T, std::size_t size>
struct is_array<Utils::Vector<T, size>> : std::true_type {};
} // namespace h5xx

namespace Writer {
namespace H5md {

using MultiArray3i = boost::multi_array<int, 3>;
using Vector1hs = Utils::Vector<hsize_t, 1>;
using Vector2hs = Utils::Vector<hsize_t, 2>;
using Vector3hs = Utils::Vector<hsize_t, 3>;
using Vector1s = Utils::Vector<size_t, 1>;
using Vector2s = Utils::Vector<size_t, 2>;
using Vector3s = Utils::Vector<size_t, 3>;

static void backup_file(const std::string &from, const std::string &to) {
  /*
   * If the file itself *and* a backup file exists, something must
   * have gone wrong.
   */
  boost::filesystem::path pfrom(from), pto(to);
  auto constexpr option_fail_if_exists = boost::filesystem::copy_options::none;
  try {
    boost::filesystem::copy_file(pfrom, pto, option_fail_if_exists);
  } catch (const boost::filesystem::filesystem_error &) {
    throw left_backupfile();
  }
}

template <typename extent_type>
static void extend_dataset(HighFive::DataSet &dataset,
                           extent_type const &change_extent) {
  auto extents = dataset.getSpace().getDimensions();
  auto const rank = extents.size();
  /* Extend the dataset for another timestep */
  for (auto i = 0u; i < rank; i++) {
    extents[i] += change_extent[i];
  }
  //std::vector<hsize_t> dst(rank);
  //std::transform(extents.begin(), extents.end(), dst.begin(),
  //		  [](long unsigned int val) {return static_cast<hsize_t>(val);});
  //H5Dset_extent(dataset.getId(), dst.data()); // extend all dims is collective
  dataset.resize(extents);
}

template <typename value_type, typename extent_type>
static void write_dataset(value_type const &data, HighFive::DataSet &dataset,
                          extent_type const &offset, extent_type const &count) {
  auto xfer_props = HighFive::DataTransferProps{};
  xfer_props.add(HighFive::UseCollectiveIO{});
  /* write the iata to the dataset. */
  dataset.select(offset, count).write(data, xfer_props);
}

template <typename value_type, typename extent_type>
static void write_dataset(value_type const &data, HighFive::DataSet &dataset,
                          extent_type const &change_extent,
                          extent_type const &offset, extent_type const &count) {
  extend_dataset(dataset, change_extent);
  /* write the data to the dataset. */
  write_dataset(data, dataset, offset, count);
}

static void write_script(std::string const &target,
                         boost::filesystem::path const &script_path) {
  if (!script_path.empty()) {
    std::ifstream scriptfile(script_path.string());
    std::string buffer((std::istreambuf_iterator<char>(scriptfile)),
                       std::istreambuf_iterator<char>());
    HighFive::File file(target, HighFive::File::Overwrite);
    file.createGroup("/parameters");
    auto group = file.createGroup("/parameters/files");
    group.createAttribute("script", buffer);
  }
}

/* Initialize the file-related variables after parameters have been set. */
void File::init_file(std::string const &file_path) {
  m_backup_filename = file_path + ".bak";
  if (m_script_path.empty()) {
    m_absolute_script_path = boost::filesystem::path();
  } else {
    boost::filesystem::path script_path(m_script_path);
    m_absolute_script_path = boost::filesystem::canonical(script_path);
  }
  auto const file_exists = boost::filesystem::exists(file_path);
  auto const backup_file_exists = boost::filesystem::exists(m_backup_filename);
  /* Perform a barrier synchronization. Otherwise one process might already
   * create the file while another still checks for its existence. */
  m_comm.barrier();
  if (file_exists) {
    if (m_h5md_specification.is_compliant(file_path)) {
      /*
       * If the file exists and has a valid H5MD structure, let's create a
       * backup of it. This has the advantage, that the new file can
       * just be deleted if the simulation crashes at some point and we
       * still have a valid trajectory backed up, from which we can restart.
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
  auto &datasets = *m_datasets;
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    if (ds.is_link)
      continue;
    datasets[ds.path()] = m_h5md_file->getDataSet(ds.path());
  }
}

void File::create_groups() {
  //HighFive::Group group = m_h5md_file->getGroup("/");
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    std::stringstream ss(ds.group);
    std::string segment;
    std::string current_path = "/";
    while (std::getline(ss, segment, '/')) {
      if (segment.empty()) continue;
      current_path += "/" + segment;
      if (!m_h5md_file->exist(current_path)) {
        m_h5md_file->createGroup(current_path);
      }
    }
  }
}

static std::vector<size_t> create_dims(hsize_t rank, hsize_t data_dim) {
  switch (rank) {
  case 3ul:
    return std::vector<size_t>{0ul, 0ul, data_dim};
  case 2ul:
    return std::vector<size_t>{0ul, data_dim};
  case 1ul:
    return std::vector<size_t>{data_dim};
  default:
    throw std::runtime_error(
        "H5MD Error: datasets with this dimension are not implemented\n");
  }
}

static std::vector<size_t> create_maxdims(hsize_t rank, hsize_t data_dim, hsize_t max_dim) {
  switch (rank) {
  case 3ul:
    return std::vector<size_t>{max_dim, max_dim, data_dim};
  case 2ul:
    return std::vector<size_t>{max_dim, max_dim};
  case 1ul:
    return std::vector<size_t>{max_dim};
  default:
    throw std::runtime_error(
        "H5MD Error: datasets with this dimension are not implemented\n");
  }
}

static std::vector<hsize_t> create_chunk_dims(hsize_t rank, hsize_t data_dim) {
  hsize_t chunk_size = (rank > 1ul) ? 1000ul : 1ul;
  switch (rank) {
  case 3ul:
    return {1ul, chunk_size, data_dim};
  case 2ul:
    return {1ul, chunk_size};
  case 1ul:
    return {chunk_size};
  default:
    throw std::runtime_error(
        "H5MD Error: datasets with this dimension are not implemented\n");
  }
}

void File::create_datasets() {
  //namespace hps = h5xx::policy::storage;
  auto &datasets = *m_datasets;
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    if (ds.is_link)
      continue;
    //auto maxdims = std::vector<hsize_t>(ds.rank, H5S_UNLIMITED);
    auto dims = create_dims(ds.rank, ds.data_dim);
    auto maxdims = create_maxdims(ds.rank, ds.data_dim, H5S_UNLIMITED);
    //auto dataspace = HighFive::DataSpace(dims, maxdims);
    auto dataspace = HighFive::DataSpace(dims, maxdims);
    //auto storage = hps::chunked(create_chunk_dims(ds.rank, ds.data_dim))
    //                   .set(hps::fill_value(-10));
    auto chunk = create_chunk_dims(ds.rank, ds.data_dim);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk));
    if (ds.type == H5T_NATIVE_INT) {
      datasets[ds.path()] =
          m_h5md_file->createDataSet<int>(ds.path(), dataspace, props); //Tmp
    } else if (ds.type == H5T_NATIVE_DOUBLE) {
      datasets[ds.path()] =
          m_h5md_file->createDataSet<double>(ds.path(), dataspace, props); //Tmp
    }
  }
}

void File::load_file(const std::string &file_path) {
  HighFive::FileAccessProps fapl;
  fapl.add(HighFive::MPIOFileAccess{m_comm, MPI_INFO_NULL});
  fapl.add(HighFive::MPIOCollectiveMetadata{});
  m_h5md_file = std::make_unique<HighFive::File>(file_path, HighFive::File::Overwrite, fapl);
  load_datasets();
}

static void write_attributes(HighFive::File &h5md_file) {
  auto h5md_group = h5md_file.createGroup("h5md");
  auto att = h5md_group.createAttribute<size_t>("version", HighFive::DataSpace::From(std::array<size_t, 2>{{1ul, 1ul}}));
  att.write(std::array<size_t, 2>{{1ul, 1ul}});
  auto h5md_creator_group = h5md_group.createGroup("creator");
  h5md_creator_group.createAttribute("name", "ESPResSo");
  h5md_creator_group.createAttribute("version", ESPRESSO_VERSION);
  auto h5md_author_group = h5md_group.createGroup("author");
  h5md_author_group.createAttribute("name", "N/A");
  auto group = h5md_file.getGroup("particles/atoms/box");
  group.createAttribute("dimension", 3);
  group.createAttribute("boundary", "periodic");
}

void File::write_units() {
  //auto const &datasets = *m_datasets;
  auto &datasets = *m_datasets;
  if (!mass_unit().empty() and (m_fields & H5MD_OUT_MASS)) {
    datasets.at("/particles/atoms/mass/value").createAttribute("unit", mass_unit());
  }
  if (!charge_unit().empty() and (m_fields & H5MD_OUT_CHARGE)) {
    datasets.at("/particles/atoms/charge/value").createAttribute("unit", charge_unit());
  }
  if (!length_unit().empty() and (m_fields & H5MD_OUT_BOX_L)) {
    datasets.at("/particles/atoms/position/value").createAttribute("unit", length_unit());
    datasets.at("/particles/atoms/box/edges/value").createAttribute("unit", length_unit());
  }
  if (!length_unit().empty() and (m_fields & H5MD_OUT_LE_OFF)) {
    datasets.at("/particles/atoms/lees_edwards/offset/value").createAttribute("unit", length_unit());
  }
  if (!velocity_unit().empty() and (m_fields & H5MD_OUT_VEL)) {
    datasets.at("/particles/atoms/velocity/value").createAttribute("unit", velocity_unit());
  }
  if (!force_unit().empty() and (m_fields & H5MD_OUT_FORCE)) {
    datasets.at("/particles/atoms/force/value").createAttribute("unit", force_unit());
  }
  if (!time_unit().empty()) {
    datasets.at("/particles/atoms/id/time").createAttribute("unit", time_unit());
  }
}

void File::create_hard_links() {
  std::string path_step = "/particles/atoms/id/step";
  std::string path_time = "/particles/atoms/id/time";
  for (auto &ds : m_h5md_specification.get_datasets()) {
    if (ds.is_link) {
      char const *from = nullptr;
      if (ds.name == "step") {
        from = path_step.c_str();
      } else if (ds.name == "time") {
        from = path_time.c_str();
      }
      assert(from != nullptr);
      //if (H5Lcreate_hard(m_h5md_file->hid(), from, m_h5md_file->hid(),
      if (H5Lcreate_hard(m_h5md_file->getId(), from, m_h5md_file->getId(),
                         ds.path().c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0) {
        throw std::runtime_error("Error creating hard link for " + ds.path());
      }
    }
  }
}

void File::create_file(const std::string &file_path) {
  if (m_comm.rank() == 0)
    write_script(file_path, m_absolute_script_path);
  m_comm.barrier();
  HighFive::FileAccessProps fapl;
  fapl.add(HighFive::MPIOFileAccess{m_comm, MPI_INFO_NULL});
  fapl.add(HighFive::MPIOCollectiveMetadata{});
  m_h5md_file = std::make_unique<HighFive::File>(file_path, HighFive::File::Overwrite, fapl);
  //m_h5md_file = std::make_unique<h5xx::file>(file_path, m_comm, MPI_INFO_NULL, h5xx::file::out);
  create_groups();
  create_datasets();
  write_attributes(*m_h5md_file);
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
    return Vector3s{1, n_part_diff, 0};
  }
  static constexpr auto count() { return Vector3s{1, 1, 3}; }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector3s{n_time_steps, prefix, 0};
  }
};

template <> struct slice_info<2> {
  static auto extent(hsize_t n_part_diff) { return Vector2s{1, n_part_diff}; }
  static constexpr auto count() { return Vector2s{1, 1}; }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector2s{n_time_steps, prefix};
  }
};

} // namespace detail

template <std::size_t dim, typename Op>
void write_td_particle_property(hsize_t prefix, hsize_t n_part_global,
                                ParticleRange const &particles,
                                HighFive::DataSet &dataset, Op op) {
  //auto const old_extents = static_cast<h5xx::dataspace>(dataset).extents();
  auto const old_extents = dataset.getSpace().getDimensions();
  auto const extent_particle_number =
      std::max(n_part_global, static_cast<hsize_t>(old_extents[1])) - old_extents[1]; //?????
  extend_dataset(dataset,
                 detail::slice_info<dim>::extent(extent_particle_number));
  auto const count = detail::slice_info<dim>::count();
  auto offset = detail::slice_info<dim>::offset(old_extents[0], prefix);
  for (auto const &p : particles) {
    //h5xx::write_dataset(dataset, op(p), h5xx::slice(offset, count));
    write_dataset(op(p), dataset, offset, count);
    // advance in the particle dimension
    offset[1] += 1;
  }
}

static void write_box(BoxGeometry const &box_geo, HighFive::DataSet &dataset) {
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1, 0});
  //h5xx::write_dataset(dataset, box_geo.length(),
  //                    h5xx::slice(Vector2hs{extents[0], 0}, Vector2hs{1, 3}));
  std::vector<size_t> offset{static_cast<size_t>(extents[0]), 0};
  std::vector<size_t> count{1ul, 3ul};
  write_dataset(box_geo.length().as_vector(), dataset, offset, count);
}

static void write_le_off(LeesEdwardsBC const &lebc, HighFive::DataSet &dataset) {
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1, 0});
  //h5xx::write_dataset(dataset, Utils::Vector<double, 1>{lebc.pos_offset},
  //                    h5xx::slice(Vector2hs{extents[0], 0}, Vector2hs{1, 1}));
  std::vector<size_t> offset{extents[0], 0};
  std::vector<size_t> count{1ul, 1ul};
  write_dataset(std::vector<double>{lebc.pos_offset}, dataset, offset, count);
}

static void write_le_dir(LeesEdwardsBC const &lebc, HighFive::DataSet &dataset) {
  auto const shear_direction = static_cast<int>(lebc.shear_direction);
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1, 0});
  //h5xx::write_dataset(dataset, Utils::Vector<int, 1>{shear_direction},
  //                    h5xx::slice(Vector2hs{extents[0], 0}, Vector2hs{1, 1}));
  std::vector<size_t> offset{extents[0], 0};
  std::vector<size_t> count{1ul, 1ul};
  write_dataset(std::vector<int>{shear_direction}, dataset, offset, count);
}

static void write_le_normal(LeesEdwardsBC const &lebc, HighFive::DataSet &dataset) {
  auto const shear_plane_normal = static_cast<int>(lebc.shear_plane_normal);
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1, 0});
  //h5xx::write_dataset(dataset, Utils::Vector<int, 1>{shear_plane_normal},
  //                    h5xx::slice(Vector2hs{extents[0], 0}, Vector2hs{1, 1}));
  std::vector<size_t> offset{extents[0], 0};
  std::vector<size_t> count{1ul, 1ul};
  write_dataset(std::vector<int>{shear_plane_normal}, dataset, offset, count);
}

void File::write(const ParticleRange &particles, double time, int step,
                 BoxGeometry const &box_geo) {
  auto &datasets = *m_datasets;
  if (m_fields & H5MD_OUT_BOX_L) {
    write_box(box_geo, datasets["/particles/atoms/box/edges/value"]);
  }
  auto const &lebc = box_geo.lees_edwards_bc();
  if (m_fields & H5MD_OUT_LE_OFF) {
    write_le_off(lebc, datasets["/particles/atoms/lees_edwards/offset/value"]);
  }
  if (m_fields & H5MD_OUT_LE_DIR) {
    write_le_dir(lebc,
                 datasets["/particles/atoms/lees_edwards/direction/value"]);
  }
  if (m_fields & H5MD_OUT_LE_NORMAL) {
    write_le_normal(lebc,
                    datasets["/particles/atoms/lees_edwards/normal/value"]);
  }

  auto const n_part_local = static_cast<int>(particles.size());
  // calculate count and offset
  int prefix = 0;
  // calculate prefix for write of the current process
  BOOST_MPI_CHECK_RESULT(MPI_Exscan,
                         (&n_part_local, &prefix, 1, MPI_INT, MPI_SUM, m_comm));

  auto const n_part_global =
      boost::mpi::all_reduce(m_comm, n_part_local, std::plus<int>());

  write_td_particle_property<2>(
      prefix, n_part_global, particles, datasets["/particles/atoms/id/value"],
      [](auto const &p) { return std::vector<int>{p.id()}; });

  {
    HighFive::DataSet &dataset = datasets["/particles/atoms/id/value"];
    auto const extents = dataset.getSpace().getDimensions();
    write_dataset(std::vector<double>{time},
                  datasets["/particles/atoms/id/time"], Vector1s{1},
                  Vector1s{extents[0]}, Vector1s{1});
    write_dataset(std::vector<int>{step},
                  datasets["/particles/atoms/id/step"], Vector1s{1},
                  Vector1s{extents[0]}, Vector1s{1});
  }

  if (m_fields & H5MD_OUT_TYPE) {
    write_td_particle_property<2>(
        prefix, n_part_global, particles,
        datasets["/particles/atoms/species/value"],
        [](auto const &p) { return std::vector<int>{p.type()}; });
  }
  if (m_fields & H5MD_OUT_MASS) {
    write_td_particle_property<2>(
        prefix, n_part_global, particles,
        datasets["/particles/atoms/mass/value"],
        [](auto const &p) { return std::vector<double>{p.mass()}; });
  }
  if (m_fields & H5MD_OUT_POS) {
    write_td_particle_property<3>(
        prefix, n_part_global, particles,
        datasets["/particles/atoms/position/value"],
        [&](auto const &p) { return box_geo.folded_position(p.pos()).as_vector(); });
  }
  if (m_fields & H5MD_OUT_IMG) {
    write_td_particle_property<3>(
        prefix, n_part_global, particles,
        datasets["/particles/atoms/image/value"], [&](auto const &p) {
          return box_geo.folded_image_box(p.pos(), p.image_box()).as_vector();
        });
  }
  if (m_fields & H5MD_OUT_VEL) {
    write_td_particle_property<3>(prefix, n_part_global, particles,
                                  datasets["/particles/atoms/velocity/value"],
                                  [](auto const &p) { return p.v().as_vector(); });
  }
  if (m_fields & H5MD_OUT_FORCE) {
    write_td_particle_property<3>(prefix, n_part_global, particles,
                                  datasets["/particles/atoms/force/value"],
                                  [](auto const &p) { return p.force().as_vector(); });
  }
  if (m_fields & H5MD_OUT_CHARGE) {
    write_td_particle_property<2>(
        prefix, n_part_global, particles,
        datasets["/particles/atoms/charge/value"],
        [](auto const &p) { return std::vector<double>{p.q()}; });
  }
  if (m_fields & H5MD_OUT_BONDS) {
    write_connectivity(particles);
  }
}

void File::write_connectivity(const ParticleRange &particles) {
  MultiArray3i bond(boost::extents[0][0][0]);
  for (auto const &p : particles) {
    auto nbonds_local = static_cast<decltype(bond)::index>(bond.shape()[1]);
    for (auto const b : p.bonds()) {
      auto const partner_ids = b.partner_ids();
      if (partner_ids.size() == 1u) {
        bond.resize(boost::extents[1][nbonds_local + 1][2]);
        bond[0][nbonds_local][0] = p.id();
        bond[0][nbonds_local][1] = partner_ids[0];
        nbonds_local++;
      }
    }
  }

  auto const n_bonds_local = static_cast<int>(bond.shape()[1]);
  auto &datasets = *m_datasets;
  int prefix_bonds = 0;
  BOOST_MPI_CHECK_RESULT(
      MPI_Exscan, (&n_bonds_local, &prefix_bonds, 1, MPI_INT, MPI_SUM, m_comm));
  auto const n_bonds_total =
      boost::mpi::all_reduce(m_comm, n_bonds_local, std::plus<int>());
  auto const extents = datasets["/connectivity/atoms/value"].getSpace().getDimensions();
  Vector3s offset_bonds = {extents[0], static_cast<size_t>(prefix_bonds), 0};
  Vector3s count_bonds = {1, static_cast<size_t>(n_bonds_local), 2};
  auto const n_bond_diff =
      std::max(static_cast<hsize_t>(n_bonds_total),
	       static_cast<hsize_t>(extents[1])) - extents[1];
  Vector3s change_extent_bonds = {1, static_cast<size_t>(n_bond_diff), 0};
  write_dataset(bond, datasets["/connectivity/atoms/value"], change_extent_bonds,
                offset_bonds, count_bonds);
}

void File::flush() { m_h5md_file->flush(); }

std::string File::file_path() const { return m_h5md_file->getName(); }

File::File(std::string file_path, std::string script_path,
           std::vector<std::string> const &output_fields, std::string mass_unit,
           std::string length_unit, std::string time_unit,
           std::string force_unit, std::string velocity_unit,
           std::string charge_unit)
    : m_script_path(std::move(script_path)), m_mass_unit(std::move(mass_unit)),
      m_length_unit(std::move(length_unit)), m_time_unit(std::move(time_unit)),
      m_force_unit(std::move(force_unit)),
      m_velocity_unit(std::move(velocity_unit)),
      m_charge_unit(std::move(charge_unit)), m_comm(boost::mpi::communicator()),
      m_fields(fields_list_to_bitfield(output_fields)),
      //m_h5md_file(std::make_unique<HighFive::File>()),
      m_datasets(std::make_unique<decltype(m_datasets)::element_type>()),
      m_h5md_specification(m_fields) {
  init_file(file_path);
}

File::~File() = default;

} /* namespace H5md */
} /* namespace Writer */
