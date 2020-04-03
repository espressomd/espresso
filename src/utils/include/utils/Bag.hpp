#ifndef ESPRESSO_UTILS_BAG_HPP
#define ESPRESSO_UTILS_BAG_HPP

#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

namespace Utils {
/**
 * @brief Bag of elements.
 *
 * A bag is a container in which the elements do not have
 * a fixed order. It can be considered an unordered variant
 * of vector, and implements the Container named requirement
 * specified in C++11.
 *
 * Elements in the container do not have a stable position and
 * removing elements can change the order of the other elements.
 * The loser contract (compared to a vector) allows removing any
 * element in the container in constant time.
 *
 * @tparam T Element Type, needs to be Swappable.
 */
template <class T> class Bag {
  /** Storage backend */
  using storage_type = std::vector<T>;

public:
  using value_type = T;
  using iterator = T *;
  using const_iterator = const T *;
  using pointer = T *;
  using reference = T &;

  /**
   * @brief Construct an empty container.
   */
  Bag() = default;

private:
  /** Underlying storage of the container */
  std::vector<T> m_storage;

  friend boost::serialization::access;
  /**
   * @brief Serialize the container.
   *
   * Serialization requires T to be serializable.
   */
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &m_storage;
  }

public:
  iterator begin() { return m_storage.data(); }
  iterator end() { return m_storage.data() + size(); }
  const_iterator begin() const { return m_storage.data(); }
  const_iterator end() const { return m_storage.data() + size(); }

  /**
   * @brief Number of elements in the container.
   */
  size_t size() const { return m_storage.size(); }

  /**
   * @brief Is the container empty?
   * @return True if there are no elements.
   */
  bool empty() const { return m_storage.empty(); }

  /**
   * @brief Capacity of the container.
   *
   * Number of elements the container can at least hold
   * without reallocating.
   */
  size_t capacity() const { return m_storage.capacity(); }

  /**
   * @brief Maximum number of elements the container can hold.
   */
  size_t max_size() const { return m_storage.max_size(); }

  /**
   * @brief Reserve storage.
   *
   * Increase capacity to at least the specified value.
   *
   *     @param new_capacity New minimum capacity.
   */
  void reserve(size_t new_capacity) { m_storage.reserve(new_capacity); }

  /**
   * @brief Resize container.
   *
   * Newly added Ts are default-initialized.
   * If the new size is larger than the capacity, all
   * iterators into the container are invalidated.
   *
   *     @param new_size Size to resize to.
   */
  void resize(size_t new_size) { m_storage.resize(new_size); }

  /**
   * @brief Remove all elements form container.
   */
  void clear() { m_storage.clear(); }

  /**
   * @brief Insert an element into the container.
   *
   * If before the call size() >= capacity(),
   * this may reallocate, in which case all
   * iterators into the container are invalidated.
   * Otherwise only the end iterator is invalidated.
   *
   * @param v Element to add.
   * @return Reference to the added element.
   */
  T &insert(T const &v) {
    m_storage.push_back(v);

    return m_storage.back();
  }
  /** @overload */
  T &insert(T &&v) {
    m_storage.push_back(std::move(v));

    return m_storage.back();
  }

  /**
   * @brief Remove element from the list.
   *
   * @param it Iterator pointing to the element to remove.
   * @return An iterator past the element that was removed.
   */
  iterator erase(iterator it) {
    *it = std::move(m_storage.back());

    m_storage.pop_back();

    return it;
  }

  /**
   * @brief Swap two Bags.
   *
   * Efficiently swap to bags by swapping
   * their contained storage.
   */
  friend void swap(Bag &lhs, Bag &rhs) {
    using std::swap;

    swap(lhs.m_storage, rhs.m_storage);
  }
};
} // namespace Utils
#endif
