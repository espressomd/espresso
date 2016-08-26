#ifndef UTILS_SERIALIZATION_ARRAY_HPP
#define UTILS_SERIALIZATION_ARRAY_HPP

#include <array>

namespace boost {
namespace serialization {
	template<typename Archive, typename T, std::size_t N>
		void serialize(Archive &ar, std::array<T, N> & a, const unsigned int) {
			ar & *static_cast<T (*)[N]>(static_cast<void *>(a.data()));
		}
}
}
#endif
