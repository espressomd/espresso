#pragma once

#define PS_FP16_INFINITY __short_as_half(0x7c00)
#define PS_FP16_NEG_INFINITY __short_as_half(0xfc00)

#define PS_FP32_INFINITY __int_as_float(0x7f800000)
#define PS_FP32_NEG_INFINITY __int_as_float(0xff800000)

#define PS_FP64_INFINITY  __longlong_as_double(0x7ff0000000000000)
#define PS_FP64_NEG_INFINITY  __longlong_as_double(0xfff0000000000000)
