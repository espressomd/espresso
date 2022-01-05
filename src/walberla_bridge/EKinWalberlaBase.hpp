#ifndef ESPRESSO_EKINWALBERLABASE_HPP
#define ESPRESSO_EKINWALBERLABASE_HPP

#include <boost/optional.hpp>
#include <string>

#include "BlockAndCell.hpp"
#include "utils/Vector.hpp"

#include "VTKHandle.hpp"
#include "blockforest/StructuredBlockForest.h"

/** Class that runs and controls the EK on WaLBerla
 */
template <typename FloatType = double> class EKinWalberlaBase {
private:
  FloatType m_diffusion;
  FloatType m_kT;
  FloatType m_valency;
  Utils::Vector3d m_ext_efield;
  bool m_advection;
  bool m_friction_coupling;

protected:
  EKinWalberlaBase(FloatType diffusion, FloatType kT, FloatType valency,
                   Utils::Vector<FloatType, 3> ext_efield, bool advection,
                   bool friction_coupling)
      : m_diffusion{diffusion}, m_kT{kT}, m_valency{valency},
        m_ext_efield{ext_efield}, m_advection{advection},
        m_friction_coupling(friction_coupling) {}

public:
  /** @brief Integrate EKin for one time step */
  virtual void integrate(const std::size_t &potential_id,
                         const std::size_t &velocity_id,
                         const std::size_t &force_id) = 0;

  /** @brief perform ghost communication of densities */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized fluxes */
  [[nodiscard]] virtual size_t stencil_size() const = 0;

  // Density
  virtual bool set_node_density(const Utils::Vector3i &node,
                                FloatType density) = 0;
  [[nodiscard]] virtual boost::optional<FloatType>
  get_node_density(const Utils::Vector3i &node) const = 0;

  [[nodiscard]] virtual bool
  set_node_flux_boundary(const Utils::Vector3i &node,
                         const Utils::Vector3d &flux) = 0;
  virtual bool remove_node_from_flux_boundary(const Utils::Vector3i &node) = 0;
  [[nodiscard]] virtual bool
  set_node_density_boundary(const Utils::Vector3i &node, double density) = 0;
  virtual bool
  remove_node_from_density_boundary(const Utils::Vector3i &node) = 0;
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_flux_boundary(const Utils::Vector3i &node,
                            bool consider_ghosts = false) const = 0;
  [[nodiscard]] virtual boost::optional<bool>
  get_node_is_density_boundary(const Utils::Vector3i &node,
                               bool consider_ghosts = false) const = 0;
  virtual void clear_flux_boundaries() = 0;
  virtual void clear_density_boundaries() = 0;

  virtual void update_flux_boundary_from_shape(std::vector<int> const &,
                                               std::vector<double> const &) = 0;
  virtual void update_flux_boundary_from_list(std::vector<int> const &,
                                              std::vector<double> const &) = 0;
  virtual void
  update_density_boundary_from_shape(std::vector<int> const &,
                                     std::vector<double> const &) = 0;
  virtual void
  update_density_boundary_from_list(std::vector<int> const &,
                                    std::vector<double> const &) = 0;

  // Global parameters
  void set_diffusion(FloatType diffusion) noexcept { m_diffusion = diffusion; }
  [[nodiscard]] FloatType get_diffusion() const noexcept { return m_diffusion; }
  void set_kT(FloatType kT) noexcept { m_kT = kT; }
  [[nodiscard]] FloatType get_kT() const noexcept { return m_kT; }
  void set_valency(FloatType valency) noexcept { m_valency = valency; }
  [[nodiscard]] FloatType get_valency() const noexcept { return m_valency; }
  void set_advection(bool advection) noexcept { m_advection = advection; }
  [[nodiscard]] bool get_advection() const noexcept { return m_advection; }
  void set_friction_coupling(bool friction_coupling) noexcept {
    m_friction_coupling = friction_coupling;
  }
  [[nodiscard]] bool get_friction_coupling() const noexcept {
    return m_friction_coupling;
  }
  [[nodiscard]] Utils::Vector<FloatType, 3> get_ext_efield() const noexcept {
    return m_ext_efield;
  }
  void set_ext_efield(const Utils::Vector<FloatType, 3> &field) noexcept {
    m_ext_efield = field;
  }

  //* @brief Fet the rng counter for thermalized LBs */
  virtual uint64_t get_rng_state() const = 0;

  /** @brief set the rng state of thermalized LBs */
  virtual void set_rng_state(uint64_t counter) = 0;

  [[nodiscard]] virtual walberla::BlockDataID get_density_id() const = 0;

  [[nodiscard]] virtual LatticeWalberla &get_lattice() const = 0;

  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *                          otherwise add a callback to write every
   *                          @p delta_N EK steps to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *                          @ref EKOutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  [[nodiscard]] virtual std::shared_ptr<VTKHandle>
  create_vtk(int delta_N, int initial_count, int flag_observables,
             std::string const &identifier, std::string const &base_folder,
             std::string const &prefix) = 0;
  /** @brief Write a VTK observable to disk.
   *
   *  @param vtk_uid          Name of the VTK object
   */
  virtual void write_vtk(std::string const &vtk_uid) = 0;
  /** @brief Toggle a VTK observable on/off.
   *
   *  @param vtk_uid          Name of the VTK object
   *  @param status           @c true to switch on, @c false to switch off
   */
  virtual void switch_vtk(std::string const &vtk_uid, bool status) = 0;

  /** @brief return a pairs of global node index and node center position */
  virtual std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts) const = 0;
  virtual ~EKinWalberlaBase() = default;
};

#endif // ESPRESSO_EKINWALBERLABASE_HPP
