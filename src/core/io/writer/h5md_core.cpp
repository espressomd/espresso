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

#include "hdf5_patches.hpp" // must appear first

#include "h5md_core.hpp"
#include "h5md_dataset.hpp"
#include "h5md_specification.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "lees_edwards/LeesEdwardsBC.hpp"

#include <config/version.hpp>

#include <utils/Vector.hpp>

#include <boost/array.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/multi_array.hpp>

#if defined(__GNUC__) or defined(__GNUG__)
// ignore false positive: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=119388
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#if __has_include(<highfive/boost.hpp>)
#include <highfive/boost.hpp>
#endif
#include <highfive/highfive.hpp>
#if defined(__GNUC__) or defined(__GNUG__)
#pragma GCC diagnostic pop
#endif

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Writer {
namespace H5md {

using MultiArray3i = boost::multi_array<int, 3>;
using Vector1hs = Utils::Vector<hsize_t, 1>;
using Vector2hs = Utils::Vector<hsize_t, 2>;
using Vector3hs = Utils::Vector<hsize_t, 3>;
using Vector1s = Utils::Vector<std::size_t, 1>;
using Vector2s = Utils::Vector<std::size_t, 2>;
using Vector3s = Utils::Vector<std::size_t, 3>;

static std::unordered_map<std::string, H5MDOutputFields> const fields_map = {
    {"all", H5MD_OUT_ALL},
    {"particle.type", H5MD_OUT_TYPE},
    {"particle.position", H5MD_OUT_POS},
    {"particle.image", H5MD_OUT_IMG},
    {"particle.velocity", H5MD_OUT_VEL},
    {"particle.force", H5MD_OUT_FORCE},
    {"particle.bonds", H5MD_OUT_BONDS},
    {"particle.charge", H5MD_OUT_CHARGE},
    {"particle.mass", H5MD_OUT_MASS},
    {"box.length", H5MD_OUT_BOX_L},
    {"lees_edwards.offset", H5MD_OUT_LE_OFF},
    {"lees_edwards.direction", H5MD_OUT_LE_DIR},
    {"lees_edwards.normal", H5MD_OUT_LE_NORMAL},
};

static auto fields_list_to_bitfield(std::vector<std::string> const &fields) {
  unsigned int bitfield = H5MD_OUT_NONE;
  for (auto const &field_name : fields) {
    if (not fields_map.contains(field_name)) {
      throw std::invalid_argument("Unknown field '" + field_name + "'");
    }
    bitfield |= fields_map.at(field_name);
  }
  return bitfield;
}

static void backup_file(std::filesystem::path const &from,
                        std::filesystem::path const &to) {
  /*
   * If the file itself *and* a backup file exists, something must
   * have gone wrong.
   */
  auto constexpr option_fail_if_exists = std::filesystem::copy_options::none;
  try {
    std::filesystem::copy_file(from, to, option_fail_if_exists);
  } catch (std::filesystem::filesystem_error const &) {
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
  dataset.resize(extents);
}

template <typename value_type, typename extent_type>
static void write_dataset(value_type const &data, HighFive::DataSet &dataset,
                          extent_type const &offset, extent_type const &count) {
  auto xfer_props = HighFive::DataTransferProps{};
  xfer_props.add(HighFive::UseCollectiveIO{});
  /* write the data to the dataset. */
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

static void write_script(HighFive::File &h5md_file,
                         std::filesystem::path const &script_path) {
  if (!script_path.empty()) {
    std::ifstream scriptfile(script_path);
    std::string buffer(std::istreambuf_iterator<char>(scriptfile),
                       std::istreambuf_iterator<char>{});
    h5md_file.createGroup("/parameters");
    auto group = h5md_file.createGroup("/parameters/files");
    group.createAttribute("script", buffer);
  }
}

/* Initialize the file-related variables after parameters have been set. */
void File::init_file() {
  auto const file_exists = std::filesystem::exists(m_file_path);
  auto const backup_file_exists = std::filesystem::exists(m_backup_path);
  /* Perform a barrier synchronization. Otherwise one process might already
   * create the file while another still checks for its existence. */
  m_comm.barrier();
  if (file_exists) {
    if (m_h5md_specification.is_compliant(m_file_path)) {
      /*
       * If the file exists and has a valid H5MD structure, let's create a
       * backup of it. This has the advantage, that the new file can
       * just be deleted if the simulation crashes at some point and we
       * still have a valid trajectory backed up, from which we can restart.
       */
      if (m_comm.rank() == 0)
        backup_file(m_file_path, m_backup_path);
      load_file();
    } else {
      throw incompatible_h5mdfile();
    }
  } else {
    if (backup_file_exists)
      throw left_backupfile();
    create_file();
  }
}

void File::load_datasets() {
  auto &datasets = m_datasets;
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    if (ds.is_link)
      continue;
    auto path = ds.path();
    datasets[path] = m_h5md_file->getDataSet(path);
  }
}

void File::create_groups() {
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    std::stringstream ss(ds.group);
    std::string segment;
    std::string current_path = "/";
    while (std::getline(ss, segment, '/')) {
      if (segment.empty())
        continue;
      current_path += "/" + segment;
      if (!m_h5md_file->exist(current_path)) {
        m_h5md_file->createGroup(current_path);
      }
    }
  }
}

static std::vector<std::size_t> create_dims(hsize_t rank, hsize_t data_dim) {
  if (rank == 3ul) {
    return {0ul, 0ul, data_dim};
  }
  if (rank == 2ul) {
    return {0ul, data_dim};
  }
  assert(rank == 1ul);
  return {0ul};
}

static std::vector<std::size_t> create_maxdims(hsize_t rank, hsize_t data_dim,
                                               hsize_t max_dim) {
  if (rank == 3ul) {
    return {max_dim, max_dim, data_dim};
  }
  if (rank == 2ul) {
    return {max_dim, max_dim};
  }
  assert(rank == 1ul);
  return {max_dim};
}

static std::vector<hsize_t> create_chunk_dims(hsize_t rank, hsize_t data_dim,
                                              hsize_t size) {
  auto const chunk_size = (rank > 1ul) ? size : hsize_t{1ul};
  if (rank == 3ul) {
    return {1ul, chunk_size, data_dim};
  }
  if (rank == 2ul) {
    return {1ul, chunk_size};
  }
  assert(rank == 1ul);
  return {chunk_size};
}

void File::create_datasets() {
  auto &datasets = m_datasets;
  for (auto const &ds : m_h5md_specification.get_datasets()) {
    if (ds.is_link)
      continue;
    auto dims = create_dims(ds.rank, ds.data_dim);
    auto maxdims = create_maxdims(ds.rank, ds.data_dim, H5S_UNLIMITED);
    auto dataspace = HighFive::DataSpace(dims, maxdims);
    auto const chunk_size = static_cast<hsize_t>(m_chunk_size);
    auto const chunk = create_chunk_dims(ds.rank, ds.data_dim, chunk_size);
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(chunk));
    auto path = ds.path();
    if (ds.type == H5T_NATIVE_INT) {
      datasets.emplace(path,
                       m_h5md_file->createDataSet<int>(path, dataspace, props));
    } else if (ds.type == H5T_NATIVE_DOUBLE) {
      datasets.emplace(
          path, m_h5md_file->createDataSet<double>(path, dataspace, props));
    }
  }
}

void File::load_file() {
  HighFive::FileAccessProps fapl;
  fapl.add(HighFive::MPIOFileAccess{m_comm, MPI_INFO_NULL});
  fapl.add(HighFive::MPIOCollectiveMetadata{});
  m_h5md_file = std::make_unique<HighFive::File>(
      m_file_path.string(), HighFive::File::ReadWrite, fapl);
  load_datasets();
}

static void write_attributes(HighFive::File &h5md_file) {
  auto h5md_group = h5md_file.createGroup("h5md");
  auto att = h5md_group.createAttribute<std::size_t>(
      "version",
      HighFive::DataSpace::From(std::array<std::size_t, 2>{{1ul, 1ul}}));
  att.write(std::array<std::size_t, 2>{{1ul, 1ul}});
  auto h5md_creator_group = h5md_group.createGroup("creator");
  h5md_creator_group.createAttribute("name", "ESPResSo");
  h5md_creator_group.createAttribute("version", ESPRESSO_VERSION);
  auto h5md_author_group = h5md_group.createGroup("author");
  h5md_author_group.createAttribute("name", "N/A");
  auto box_path = "/particles/atoms/box";
  if (h5md_file.exist(box_path)) {
    auto group = h5md_file.getGroup(box_path);
    group.createAttribute("dimension", 3);
    group.createAttribute("boundary", "periodic");
  }
}

void File::write_units() {
  auto &datasets = m_datasets;
  if (!mass_unit().empty() and (m_fields & H5MD_OUT_MASS)) {
    datasets.at("/particles/atoms/mass/value")
        .createAttribute("unit", mass_unit());
  }
  if (!charge_unit().empty() and (m_fields & H5MD_OUT_CHARGE)) {
    datasets.at("/particles/atoms/charge/value")
        .createAttribute("unit", charge_unit());
  }
  if (!length_unit().empty() and (m_fields & H5MD_OUT_BOX_L)) {
    datasets.at("/particles/atoms/position/value")
        .createAttribute("unit", length_unit());
    datasets.at("/particles/atoms/box/edges/value")
        .createAttribute("unit", length_unit());
  }
  if (!length_unit().empty() and (m_fields & H5MD_OUT_LE_OFF)) {
    datasets.at("/particles/atoms/lees_edwards/offset/value")
        .createAttribute("unit", length_unit());
  }
  if (!velocity_unit().empty() and (m_fields & H5MD_OUT_VEL)) {
    datasets.at("/particles/atoms/velocity/value")
        .createAttribute("unit", velocity_unit());
  }
  if (!force_unit().empty() and (m_fields & H5MD_OUT_FORCE)) {
    datasets.at("/particles/atoms/force/value")
        .createAttribute("unit", force_unit());
  }
  if (!time_unit().empty()) {
    datasets.at("/particles/atoms/id/time")
        .createAttribute("unit", time_unit());
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
      if (H5Lcreate_hard(m_h5md_file->getId(), from, m_h5md_file->getId(),
                         ds.path().c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0) {
        throw std::runtime_error("Error creating hard link for " + ds.path());
      }
    }
  }
}

void File::create_file() {
  HighFive::FileAccessProps fapl;
  fapl.add(HighFive::MPIOFileAccess{m_comm, MPI_INFO_NULL});
  fapl.add(HighFive::MPIOCollectiveMetadata{});
  m_h5md_file = std::make_unique<HighFive::File>(m_file_path.string(),
                                                 HighFive::File::Create, fapl);
  write_script(*m_h5md_file, m_absolute_script_path);
  create_groups();
  create_datasets();
  write_attributes(*m_h5md_file);
  write_units();
  create_hard_links();
}

void File::close() {
  // release datasets first, since they reference the open file
  m_datasets.clear();
  // close the file in parallel (this also flushes the superblock)
  m_h5md_file.reset();
  // wait for all ranks to be successful before deleting the backup file
  m_comm.barrier();
  if (m_comm.rank() == 0) {
    std::filesystem::remove(m_backup_path);
  }
}

namespace detail {

template <std::size_t rank> struct slice_info {};

template <> struct slice_info<3> {
  static auto extent(hsize_t n_part_diff) {
    return Vector3s{1ul, n_part_diff, 0ul};
  }
  static constexpr auto count(std::size_t local_n_part) {
    return Vector3s{1ul, local_n_part, 3ul};
  }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector3s{n_time_steps, prefix, 0ul};
  }
  template <typename T>
  static boost::multi_array<T, 3> reshape(std::vector<T> const &v1d,
                                          Vector3s const &count) {
    if (v1d.empty()) {
      boost::multi_array<T, 3> data(boost::extents[0][0][0]);
      return data;
    }
    auto const rows = count[1];
    auto const cols = count[2];

    boost::multi_array<T, 3> data(
        boost::extents[1][static_cast<long>(rows)][static_cast<long>(cols)]);

    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        data[0][i][j] = v1d[cols * i + j];
      }
    }

    return data;
  }
};

template <> struct slice_info<2> {
  static auto extent(hsize_t n_part_diff) { return Vector2s{1ul, n_part_diff}; }
  static constexpr auto count(std::size_t local_n) {
    return Vector2s{1ul, local_n};
  }
  static auto offset(hsize_t n_time_steps, hsize_t prefix) {
    return Vector2s{n_time_steps, prefix};
  }
  template <typename T>
  static boost::multi_array<T, 2> reshape(std::vector<T> const &v1d,
                                          Vector2s const &count) {
    if (v1d.empty()) {
      boost::multi_array<T, 2> data(boost::extents[0][0]);
      return data;
    }
    auto const cols = count[1];

    boost::multi_array<T, 2> data(boost::extents[1][static_cast<long>(cols)]);

    for (std::size_t i = 0; i < cols; ++i) {
      data[0][i] = v1d[i];
    }

    return data;
  }
};

template <typename T> struct get_buffer_traits {};

template <typename T>
  requires std::is_arithmetic_v<T>
struct get_buffer_traits<T> {
  using type = T;
  constexpr static std::size_t dim = 1ul;
};

template <typename T, std::size_t N>
  requires std::is_arithmetic_v<T>
struct get_buffer_traits<Utils::Vector<T, N>> {
  using type = T;
  constexpr static std::size_t dim = N;
};

template <typename Functor> class ParticleDataSerializer {
  using RetVal = std::decay_t<std::invoke_result_t<Functor, Particle const &>>;
  Functor m_getter;

  template <typename T>
    requires std::is_arithmetic_v<T>
  void serialize(auto &buffer, T const &value) const {
    buffer.emplace_back(value);
  }

  template <typename T, std::size_t N>
  void serialize(auto &buffer, Utils::Vector<T, N> const &value) const {
    buffer.insert(buffer.end(), value.cbegin(), value.cend());
  }

public:
  explicit ParticleDataSerializer(Functor lambda) : m_getter{lambda} {}

  auto operator()(ParticleRange const &particles) const {
    auto constexpr value_dim = get_buffer_traits<RetVal>::dim;
    std::vector<typename get_buffer_traits<RetVal>::type> buffer{};
    buffer.reserve(particles.size() * value_dim);
    for (auto const &p : particles) {
      serialize(buffer, m_getter(p));
    }
    return buffer;
  }
};

template <typename Functor> auto make_serializer(Functor lambda) {
  return ParticleDataSerializer<Functor>{lambda};
}
template <typename RetVal>
auto make_serializer(RetVal (Particle::*getter)() const) {
  auto kernel = [getter](Particle const &p) -> RetVal { return (p.*getter)(); };
  return ParticleDataSerializer<decltype(kernel)>{std::move(kernel)};
}

} // namespace detail

template <std::size_t dim, typename Serializer>
void write_td_particle_property(hsize_t prefix, hsize_t n_part_global,
                                ParticleRange const &particles,
                                HighFive::DataSet &dataset,
                                Serializer serializer) {
  auto const n_part_local = static_cast<hsize_t>(particles.size());
  auto const old_extents = dataset.getSpace().getDimensions();
  auto const extent_n_part =
      std::max(n_part_global, static_cast<hsize_t>(old_extents[1])) -
      old_extents[1];
  extend_dataset(dataset, detail::slice_info<dim>::extent(extent_n_part));
  auto const count = detail::slice_info<dim>::count(n_part_local);
  auto const offset = detail::slice_info<dim>::offset(old_extents[0], prefix);
  HighFive::DataType dtype = dataset.getDataType();
  auto buffer = serializer(particles);
  write_dataset(detail::slice_info<dim>::reshape(buffer, count), dataset,
                offset, count);
}

static void write_box(BoxGeometry const &box_geo, HighFive::DataSet &dataset) {
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1ul, 0ul});
  Vector2s const offset{extents[0], 0ul};
  Vector2s const count{1ul, 3ul};
  auto const data = box_geo.length().as_vector();
  write_dataset(detail::slice_info<2>::reshape(data, count), dataset, offset,
                count);
}

static void write_le_off(LeesEdwardsBC const &lebc,
                         HighFive::DataSet &dataset) {
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1ul, 0ul});
  Vector2s const offset{extents[0], 0ul};
  Vector2s const count{1ul, 1ul};
  auto const data = std::vector<double>{lebc.pos_offset};
  write_dataset(detail::slice_info<2>::reshape(data, count), dataset, offset,
                count);
}

static void write_le_dir(LeesEdwardsBC const &lebc,
                         HighFive::DataSet &dataset) {
  auto const shear_direction = static_cast<int>(lebc.shear_direction);
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1ul, 0ul});
  Vector2s const offset{extents[0], 0ul};
  Vector2s const count{1ul, 1ul};
  auto const data = std::vector<int>{shear_direction};
  write_dataset(detail::slice_info<2>::reshape(data, count), dataset, offset,
                count);
}

static void write_le_normal(LeesEdwardsBC const &lebc,
                            HighFive::DataSet &dataset) {
  auto const shear_plane_normal = static_cast<int>(lebc.shear_plane_normal);
  auto const extents = dataset.getSpace().getDimensions();
  extend_dataset(dataset, Vector2hs{1ul, 0ul});
  Vector2s const offset{extents[0], 0ul};
  Vector2s const count{1ul, 1ul};
  auto const data = std::vector<int>{shear_plane_normal};
  write_dataset(detail::slice_info<2>::reshape(data, count), dataset, offset,
                count);
}

void File::write(const ParticleRange &particles, double time, int step,
                 BoxGeometry const &box_geo) {
  auto &datasets = m_datasets;
  if (m_fields & H5MD_OUT_BOX_L) {
    write_box(box_geo, datasets.at("/particles/atoms/box/edges/value"));
  }
  auto const &lebc = box_geo.lees_edwards_bc();
  if (m_fields & H5MD_OUT_LE_OFF) {
    write_le_off(lebc,
                 datasets.at("/particles/atoms/lees_edwards/offset/value"));
  }
  if (m_fields & H5MD_OUT_LE_DIR) {
    write_le_dir(lebc,
                 datasets.at("/particles/atoms/lees_edwards/direction/value"));
  }
  if (m_fields & H5MD_OUT_LE_NORMAL) {
    write_le_normal(lebc,
                    datasets.at("/particles/atoms/lees_edwards/normal/value"));
  }

  // calculate particle count and offset
  static_assert(sizeof(hsize_t) == 8ul);
  auto const n_part_local = static_cast<hsize_t>(particles.size());
  hsize_t prefix{0ul};
  BOOST_MPI_CHECK_RESULT(
      MPI_Exscan, (&n_part_local, &prefix, 1, MPI_UINT64_T, MPI_SUM, m_comm));
  auto const n_part_global =
      boost::mpi::all_reduce(m_comm, n_part_local, std::plus<hsize_t>());

  write_td_particle_property<2>(prefix, n_part_global, particles,
                                datasets.at("/particles/atoms/id/value"),
                                detail::make_serializer(&Particle::id));

  {
    auto const time_ext =
        datasets.at("/particles/atoms/id/time").getSpace().getDimensions()[0];
    write_dataset(std::vector<double>{time},
                  datasets.at("/particles/atoms/id/time"), Vector1s{1},
                  Vector1s{time_ext}, Vector1s{1});
    auto const step_ext =
        datasets.at("/particles/atoms/id/step").getSpace().getDimensions()[0];
    write_dataset(std::vector<int>{step},
                  datasets.at("/particles/atoms/id/step"), Vector1s{1},
                  Vector1s{step_ext}, Vector1s{1});
  }

  if (m_fields & H5MD_OUT_TYPE) {
    write_td_particle_property<2>(prefix, n_part_global, particles,
                                  datasets.at("/particles/atoms/species/value"),
                                  detail::make_serializer(&Particle::type));
  }
  if (m_fields & H5MD_OUT_MASS) {
    write_td_particle_property<2>(prefix, n_part_global, particles,
                                  datasets.at("/particles/atoms/mass/value"),
                                  detail::make_serializer(&Particle::mass));
  }
  if (m_fields & H5MD_OUT_POS) {
    write_td_particle_property<3>(
        prefix, n_part_global, particles,
        datasets.at("/particles/atoms/position/value"),
        detail::make_serializer([&](Particle const &p) {
          return box_geo.folded_position(p.pos());
        }));
  }
  if (m_fields & H5MD_OUT_IMG) {
    write_td_particle_property<3>(
        prefix, n_part_global, particles,
        datasets.at("/particles/atoms/image/value"),
        detail::make_serializer([&](Particle const &p) {
          return box_geo.folded_image_box(p.pos(), p.image_box());
        }));
  }
  if (m_fields & H5MD_OUT_VEL) {
    write_td_particle_property<3>(
        prefix, n_part_global, particles,
        datasets.at("/particles/atoms/velocity/value"),
        detail::make_serializer(&Particle::v));
  }
  if (m_fields & H5MD_OUT_FORCE) {
    write_td_particle_property<3>(prefix, n_part_global, particles,
                                  datasets.at("/particles/atoms/force/value"),
                                  detail::make_serializer(&Particle::force));
  }
  if (m_fields & H5MD_OUT_CHARGE) {
    write_td_particle_property<2>(prefix, n_part_global, particles,
                                  datasets.at("/particles/atoms/charge/value"),
                                  detail::make_serializer(&Particle::q));
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
      auto const &partner_ids = b.partner_ids();
      if (partner_ids.size() == 1u) {
        bond.resize(boost::extents[1][nbonds_local + 1][2]);
        bond[0][nbonds_local][0] = p.id();
        bond[0][nbonds_local][1] = partner_ids[0];
        nbonds_local++;
      }
    }
  }

  auto const n_bonds_local = static_cast<int>(bond.shape()[1]);
  auto &datasets = m_datasets;
  int prefix_bonds = 0;
  BOOST_MPI_CHECK_RESULT(
      MPI_Exscan, (&n_bonds_local, &prefix_bonds, 1, MPI_INT, MPI_SUM, m_comm));
  auto const n_bonds_total =
      boost::mpi::all_reduce(m_comm, n_bonds_local, std::plus<int>());
  auto const extents =
      datasets.at("/connectivity/atoms/value").getSpace().getDimensions();
  Vector3s offset_bonds = {extents[0], static_cast<std::size_t>(prefix_bonds),
                           0};
  Vector3s count_bonds = {1, static_cast<std::size_t>(n_bonds_local), 2};
  auto const n_bond_diff = std::max(static_cast<hsize_t>(n_bonds_total),
                                    static_cast<hsize_t>(extents[1])) -
                           extents[1];
  Vector3s change_extent_bonds = {1, static_cast<std::size_t>(n_bond_diff), 0};
  write_dataset(bond, datasets.at("/connectivity/atoms/value"),
                change_extent_bonds, offset_bonds, count_bonds);
}

void File::flush() { m_h5md_file->flush(); }

std::vector<std::string> File::valid_fields() const {
  auto const view = std::views::elements<0>(fields_map);
  return {view.begin(), view.end()};
}

File::File(std::filesystem::path file_path, std::filesystem::path script_path,
           std::vector<std::string> const &output_fields, std::string mass_unit,
           std::string length_unit, std::string time_unit,
           std::string force_unit, std::string velocity_unit,
           std::string charge_unit, int chunk_size)
    : m_file_path(std::move(file_path)), m_backup_path(m_file_path),
      m_script_path(std::move(script_path)),
      m_absolute_script_path(
          m_script_path.empty()
              ? std::filesystem::path()
              : std::filesystem::weakly_canonical(m_script_path)),
      m_mass_unit(std::move(mass_unit)), m_length_unit(std::move(length_unit)),
      m_time_unit(std::move(time_unit)), m_force_unit(std::move(force_unit)),
      m_velocity_unit(std::move(velocity_unit)),
      m_charge_unit(std::move(charge_unit)), m_chunk_size(chunk_size),
      m_comm(boost::mpi::communicator()),
      m_fields(fields_list_to_bitfield(output_fields)), m_datasets(),
      m_h5md_specification(m_fields) {
  if (chunk_size <= 0) {
    throw std::domain_error("Parameter 'chunk_size' must be > 0");
  }
  m_backup_path += ".bak";
  init_file();
}

File::~File() {
  if (m_h5md_file) {
    close();
  }
}

} /* namespace H5md */
} /* namespace Writer */
