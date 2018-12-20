#ifndef UTILS_INTEGER_SEQUENCE
#define UTILS_INTEGER_SEQUENCE

namespace Utils {

template<class T, T... I>
struct integer_sequence {
  using value_type = T;
  static constexpr size_t size() noexcept { return sizeof...(I); }
};

template<size_t... I>
using index_sequence = integer_sequence<size_t, I...>;

namespace detail {

template <typename Seq, size_t SeqSize, size_t Rem>
struct Extend;

template <typename T, T... Ints, size_t SeqSize>
struct Extend<integer_sequence<T, Ints...>, SeqSize, 0> {
  using type = integer_sequence<T, Ints..., (Ints + SeqSize)...>;
};

template <typename T, T... Ints, size_t SeqSize>
struct Extend<integer_sequence<T, Ints...>, SeqSize, 1> {
  using type = integer_sequence<T, Ints..., (Ints + SeqSize)..., 2 * SeqSize>;
};

template <typename T, size_t N>
struct Gen {
  using type =
      typename Extend<typename Gen<T, N / 2>::type, N / 2, N % 2>::type;
};

template <typename T>
struct Gen<T, 0> {
  using type = integer_sequence<T>;
};

}

template <typename T, T N>
using make_integer_sequence = typename detail::Gen<T, N>::type;

template <size_t N>
using make_index_sequence = make_integer_sequence<size_t, N>;

template <typename... Ts>
using index_sequence_for = make_index_sequence<sizeof...(Ts)>;

}
#endif
