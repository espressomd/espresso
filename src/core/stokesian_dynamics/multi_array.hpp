#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <type_traits>

#ifdef __CUDACC__
#define DEVICE_FUNC __host__ __device__
#else
#define DEVICE_FUNC
#endif

#if __cplusplus == 201402L
#define CXX14_CONSTEXPR constexpr
#else
#define CXX14_CONSTEXPR
#endif

/// \cond
namespace meta {

// product

template <typename T, T...>
struct product;

template <typename T, T dim0, T... dim>
struct product<T, dim0, dim...> {
    static constexpr T const value = dim0 * product<T, dim...>::value;
};

template <typename T>
struct product<T> {
    static constexpr T const value = 1;
};

// all

template <bool...>
struct all;

template <bool... b>
struct all<true, b...> {
    static constexpr bool const value = true && all<b...>::value;
};

template <bool... b>
struct all<false, b...> {
    static constexpr bool const value = false && all<b...>::value;
};

template <>
struct all<> {
    static constexpr bool const value = true;
};

// linearized_index

template <size_t n, size_t... N>
struct linearized_index {
    template <typename... Idx>
    DEVICE_FUNC constexpr std::size_t operator()(Idx... idx) const noexcept {
        using unpack = std::size_t[];
        return unpack{std::size_t(idx)...}[n] +
               unpack{std::size_t(N)...}[n] *
                   linearized_index<n - 1, N...>{}(idx...);
    }
};

template <size_t... N>
struct linearized_index<0, N...> {
    template <typename... Idx>
    DEVICE_FUNC constexpr std::size_t operator()(Idx... idx) const noexcept {
        using unpack = std::size_t[];
        return unpack{std::size_t(idx)...}[0];
    }
};

// check_bounds

template <size_t n, size_t... N>
struct check_bounds {
    template <typename... Idx>
    constexpr bool operator()(Idx... idx) const {
        using unpack = std::size_t[];
        return unpack{std::size_t(idx)...}[n] < unpack{std::size_t(N)...}[n]
                   ? check_bounds<n - 1, N...>{}(idx...)
                   : throw std::out_of_range("index out of bounds: " +
                                             std::to_string(n));
    }
};

template <size_t... N>
struct check_bounds<0, N...> {
    template <typename... Idx>
    constexpr bool operator()(Idx... idx) const {
        using unpack = std::size_t[];
        return unpack{std::size_t(idx)...}[0] < unpack{std::size_t(N)...}[0]
                   ? false
                   : throw std::out_of_range("index out of bounds: " +
                                             std::to_string(0));
    }
};

} // namespace meta
/// \endcond

/// Multidimensional array similar to `std::array`
///
/// Storage is always in row-major order for easy list initialization
///
/// \tparam T value type
/// \tparam N list of dimensions
template <typename T, std::size_t... N>
class multi_array {
public:
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = value_type &;
    using const_reference = const value_type &;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using iterator = pointer;

private:
    template <size_type... s>
    using size_product = meta::product<size_type, s...>;

    // Storage
    value_type m_data[size_product<N...>::value];

public:
    DEVICE_FUNC constexpr multi_array() {}

    /// Constructor for aggregate initializers
    template <typename... U>
    DEVICE_FUNC constexpr multi_array(U... data)
        : m_data{value_type(data)...} {}

    /// fill the container with specified value
    DEVICE_FUNC void fill(value_type const &value) noexcept(
        std::is_nothrow_copy_assignable<value_type>::value) {
        for (iterator it = begin(); it != end(); ++it) {
            *it = value;
        }
    }

    /// swaps the contents
    DEVICE_FUNC void swap(multi_array &other) noexcept(
        std::is_nothrow_move_assignable<value_type>::value) {
        iterator oit = other.begin();
        for (iterator it = begin(); it != end(); ++it, (void)++oit) {
            value_type tmp = std::move(*it);
            *it = std::move(*oit);
            *oit = std::move(tmp);
        }
    }

    /// returns an iterator to the beginning
    DEVICE_FUNC constexpr iterator begin() const noexcept {
        return iterator(data());
    }

    /// returns an iterator to the end
    DEVICE_FUNC constexpr iterator end() const noexcept {
        return iterator(data() + size_product<N...>::value);
    }

    /// direct access to the underlying array
    DEVICE_FUNC CXX14_CONSTEXPR pointer data() noexcept { return m_data; }
    /// \overload DEVICE_FUNC CXX14_CONSTEXPR pointer data()
    DEVICE_FUNC CXX14_CONSTEXPR const_pointer data() const noexcept {
        return m_data;
    }

    /// access specified element
    ///
    /// \throws std::out_of_range if index is out of range (only in DEBUG mode)
    template <typename... Idx>
    DEVICE_FUNC CXX14_CONSTEXPR reference operator()(Idx... idx)
#ifdef NDEBUG
        noexcept
#endif
    {
        static_assert(sizeof...(idx) == sizeof...(N), "dimension mismatch");
        static_assert(
            meta::all<std::is_convertible<Idx, size_type>::value...>::value,
            "type mismatch");
        return
#if !defined(NDEBUG) && !defined(__CUDACC__)
            meta::check_bounds<sizeof...(idx) - 1, N...>{}(idx...),
#endif
            m_data[meta::linearized_index<sizeof...(idx) - 1, N...>{}(idx...)];
    }

    /// \overload DEVICE_FUNC CXX14_CONSTEXPR reference operator()(Idx... idx)
    template <typename... Idx>
    DEVICE_FUNC constexpr const_reference operator()(Idx... idx) const
#ifdef NDEBUG
        noexcept
#endif
    {
        static_assert(sizeof...(idx) == sizeof...(N), "dimension mismatch");
        static_assert(
            meta::all<std::is_convertible<Idx, size_type>::value...>::value,
            "type mismatch");
        return
#if !defined(NDEBUG) && !defined(__CUDACC__)
            meta::check_bounds<sizeof...(idx) - 1, N...>{}(idx...),
#endif
            m_data[meta::linearized_index<sizeof...(idx) - 1, N...>{}(idx...)];
    }

    /// returns the size along dimension \p i
    ///
    /// \tparam i dimension
    template <size_type i>
    DEVICE_FUNC constexpr size_type size() const noexcept {
        static_assert(i < sizeof...(N), "index out of bounds");
        using unpack = std::size_t[];
        return unpack{std::size_t(N)...}[i];
    }

    /// returns the number of elements
    DEVICE_FUNC constexpr size_type size() const noexcept {
        return size_product<N...>::value;
    }
};

template <typename T, std::size_t M, std::size_t N>
DEVICE_FUNC CXX14_CONSTEXPR multi_array<T,M,N> outer(multi_array<T,M> const & a, multi_array<T,N> const & b) {
    multi_array<T,M,N> c;
    for (std::size_t i = 0; i < M; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            c(i, j) = a(i) * b(j);
        }
    }
    return c;
}

#undef CXX14_CONSTEXPR
#undef DEVICE_FUNC
