#ifndef OBJECT_REGISTRY_HPP
#define OBJECT_REGISTRY_HPP

#include <algorithm>

template <typename Container> class ObjectRegistry {
public:
  typedef Container container_type;
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;
  using size_type = typename Container::size_type;

  void add(typename container_type::value_type const &c) {
    m_container.push_back(c);
  }

  void remove(typename container_type::value_type const &c) {
    auto it = std::find(std::begin(m_container), std::end(m_container), c);

    if (it != std::end(m_container)) {
      m_container.erase(it);
    }
  }

  bool empty() const { return m_container.empty(); }
  size_type size() const { return m_container.size(); }

  const_iterator begin() const { return std::begin(m_container); }
  const_iterator end() const { return std::end(m_container); }

private:
  Container m_container;
};

#endif
