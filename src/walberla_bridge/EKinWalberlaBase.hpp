#ifndef ESPRESSO_EKINWALBERLABASE_HPP
#define ESPRESSO_EKINWALBERLABASE_HPP

#include <boost/optional.hpp>
#include <string>

#include "utils/Vector.hpp"

#include "blockforest/StructuredBlockForest.h"

#include "WalberlaBlockForest.hpp"

/** Class that runs and controls the LB on WaLBerla
 */
template <typename FloatType = double> class EKinWalberlaBase {
private:
  FloatType m_diffusion;
  FloatType m_kT;

protected:
  EKinWalberlaBase(FloatType diffusion, FloatType kT)
      : m_diffusion{diffusion}, m_kT{kT} {}

public:
  /** @brief Integrate EKin for one time step */
  virtual void integrate() = 0;

  /** @brief perform ghost communication of densities */
  virtual void ghost_communication() = 0;

  /** @brief Number of discretized fluxes */
  [[nodiscard]] virtual size_t stencil_size() const = 0;

  // Density
  virtual bool set_node_density(const Utils::Vector3i &node,
                                FloatType density) = 0;
  [[nodiscard]] virtual boost::optional<FloatType>
  get_node_density(const Utils::Vector3i &node) const = 0;

  //  virtual bool remove_node_from_boundary(const Utils::Vector3i &node) = 0;
  //  [[nodiscard]] virtual boost::optional<bool>
  //  get_node_is_boundary(const Utils::Vector3i &node,
  //                       bool consider_ghosts = false) const = 0;
  //  virtual void clear_boundaries() = 0;

  // Global parameters
  void set_diffusion(double diffusion) { m_diffusion = diffusion; }
  [[nodiscard]] FloatType get_diffusion() const { return m_diffusion; }
  void set_kT(double kT) { m_kT = kT; }
  [[nodiscard]] FloatType get_kT() const { return m_kT; }

  //* @brief Fet the rng counter for thermalized LBs */
  virtual uint64_t get_rng_state() const = 0;

  /** @brief set the rng state of thermalized LBs */
  virtual void set_rng_state(uint64_t counter) = 0;

  [[nodiscard]] virtual const walberla::WalberlaBlockForest *
  get_blockforest() const = 0;
  /** @brief Create a VTK observable.
   *
   *  @param delta_N          Write frequency, if 0 write a single frame,
   *         otherwise add a callback to write every @p delta_N LB steps
   *         to a new file
   *  @param initial_count    Initial execution count
   *  @param flag_observables Which observables to measure (OR'ing of
   *         @ref OutputVTK values)
   *  @param identifier       Name of the VTK dataset
   *  @param base_folder      Path to the VTK folder
   *  @param prefix           Prefix of the VTK files
   */
  virtual void create_vtk(unsigned delta_N, unsigned initial_count,
                          unsigned flag_observables,
                          std::string const &identifier,
                          std::string const &base_folder,
                          std::string const &prefix) = 0;
  /** @brief Write a VTK observable to disk.
   *
   *  @param vtk_uid          Name of the VTK object
   */
  virtual void write_vtk(std::string const &vtk_uid) = 0;
  /** @brief Toggle a VTK observable on/off.
   *
   *  @param vtk_uid          Name of the VTK object
   *  @param status           1 to switch on, 0 to switch off
   */
  virtual void switch_vtk(std::string const &vtk_uid, int status) = 0;

  /** @brief return a pairs of global node index and node center position */
  virtual std::vector<std::pair<Utils::Vector3i, Utils::Vector3d>>
  node_indices_positions(bool include_ghosts) const = 0;
  virtual ~EKinWalberlaBase() = default;
};

/** @brief EKin statistics to write to VTK files */
enum class OutputVTK : unsigned {
  density = 1u << 0u,
};

#endif // ESPRESSO_EKINWALBERLABASE_HPP
