#ifndef UTILS_ATTRIBUTES_HPP
#define UTILS_ATTRIBUTES_HPP

/* __GNUC__ catches gcc, intel and clang,
 * MSVS does not have the attribute */
#if defined(__GNUC__)
#define FUNCTION_ATTRIBUTE_CONST __attribute__((const))
#else
#define FUNCTION_ATTRIBUTE_CONST
#endif
#endif
