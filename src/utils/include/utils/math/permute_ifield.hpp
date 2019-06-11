#ifndef UTILS_MATH_PERMUTE_IFIELD_HPP
#define UTILS_MATH_PERMUTE_IFIELD_HPP

namespace Utils {
/** permute an integer array field of size size about permute positions. */
inline void permute_ifield(int *field, int size, int permute) {
  if (permute == 0)
    return;
  if (permute < 0)
    permute = (size + permute);
  while (permute > 0) {
    int tmp = field[0];
    for (int i = 1; i < size; i++)
      field[i - 1] = field[i];
    field[size - 1] = tmp;
    permute--;
  }
}
} // namespace Utils

#endif
