#ifndef ESPRESSO_EKCONTAINER_HPP
#define ESPRESSO_EKCONTAINER_HPP

#ifdef LB_WALBERLA

#include "electrokinetics/PoissonSolver/PoissonSolver.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

template <class EKSpecies> class EKContainer {
  using container_type = std::vector<std::shared_ptr<EKSpecies>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  container_type m_ekcontainer;
  double m_tau{};

  // TODO: this could be moved to the Scriptinterface
  std::shared_ptr<walberla::PoissonSolver> m_poissonsolver;

public:
  void add(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) ==
           m_ekcontainer.end());

    m_ekcontainer.emplace_back(c);
  }
  void remove(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) !=
           m_ekcontainer.end());
    m_ekcontainer.erase(
        std::remove(m_ekcontainer.begin(), m_ekcontainer.end(), c),
        m_ekcontainer.end());
  }

  iterator begin() noexcept { return m_ekcontainer.begin(); }
  iterator end() noexcept { return m_ekcontainer.end(); }
  const_iterator begin() const noexcept { return m_ekcontainer.begin(); }
  const_iterator end() const noexcept { return m_ekcontainer.end(); }
  [[nodiscard]] bool empty() const noexcept { return m_ekcontainer.empty(); }

  void on_boxl_change() const {
    if (not this->empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active EKs.");
    }
  }

  void set_poissonsolver(
      std::shared_ptr<walberla::PoissonSolver> const &solver) noexcept {
    m_poissonsolver = solver;
  }

  [[nodiscard]] bool is_poissonsolver_set() const noexcept {
    return m_poissonsolver != nullptr;
  }

  [[nodiscard]] double get_tau() const noexcept { return m_tau; }
  void set_tau(double tau) noexcept { m_tau = tau; }

  void reset_charge() const { m_poissonsolver->reset_charge_field(); }
  void add_charge(const walberla::domain_decomposition::BlockDataID &id,
                  double valency, bool is_double_precision) const {
    m_poissonsolver->add_charge_to_field(id, valency, is_double_precision);
  }

  void solve_poisson() const { m_poissonsolver->solve(); }

  [[nodiscard]] walberla::domain_decomposition::BlockDataID
  get_potential_field_id() const {
    return m_poissonsolver->get_potential_field_id();
  }
};

#endif // LB_WALBERLA

#endif // ESPRESSO_EKCONTAINER_HPP