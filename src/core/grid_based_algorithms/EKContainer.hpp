#ifndef ESPRESSO_EKCONTAINER_HPP
#define ESPRESSO_EKCONTAINER_HPP

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

public:
  void add(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) ==
           m_ekcontainer.end());

    m_ekcontainer.emplace_back(c);
    // TODO: callback on_ek_change?
  }
  void remove(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) !=
           m_ekcontainer.end());
    m_ekcontainer.erase(
        std::remove(m_ekcontainer.begin(), m_ekcontainer.end(), c),
        m_ekcontainer.end());
    // TODO: callback on_ek_change?
  }

  iterator begin() { return m_ekcontainer.begin(); }
  iterator end() { return m_ekcontainer.end(); }
  const_iterator begin() const { return m_ekcontainer.begin(); }
  const_iterator end() const { return m_ekcontainer.end(); }

  void on_boxl_change() const {
    if (not this->empty()) {
      throw std::runtime_error("The box size can not be changed because there "
                               "are active EKs.");
    }
  }
};

#endif // ESPRESSO_EKCONTAINER_HPP