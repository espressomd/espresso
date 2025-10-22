/*
 * Copyright (C) 2024 The ESPResSo project
 * Copyright (C) 2024 The waLBerla project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* section */

#if defined(__NVCC__)
{% for name in defines -%}
#define {{name}} {{defines_table["nvcc"][name]}}
{% endfor -%}
{% if pragmas["nvcc"] -%}
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
{% for pragma in pragmas["nvcc"] -%}
#pragma nv_{{pragma}}
{% endfor -%}
#else
#pragma push
{% for pragma in pragmas["nvcc"] -%}
#pragma {{pragma}}
{% endfor -%}
#endif // defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
{% endif -%}
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
{% for name in defines -%}
#define {{name}} {{defines_table["clang_dev"][name]}}
{% endfor -%}
{% if pragmas["clang_dev"] -%}
#pragma clang diagnostic push
{% for pragma in pragmas["clang_dev"] -%}
#pragma clang diagnostic ignored "{{pragma}}"
{% endfor -%}
{% endif -%}
#else
// clang compiling CUDA code in host mode
{% for name in defines -%}
#define {{name}} {{defines_table["clang_host"][name]}}
{% endfor -%}
{% if pragmas["clang_host"] -%}
#pragma clang diagnostic push
{% for pragma in pragmas["clang_host"] -%}
#pragma clang diagnostic ignored "{{pragma}}"
{% endfor -%}
{% endif -%}
#endif // defined(__CUDA_ARCH__)
#endif // defined(__CUDA__)
#elif defined(__GNUC__) or defined(__GNUG__)
{% for name in defines -%}
#define {{name}} {{defines_table["gcc"][name]}}
{% endfor -%}
{% if pragmas["gcc"] -%}
#pragma GCC diagnostic push
{% for pragma in pragmas["gcc"] -%}
#pragma GCC diagnostic ignored "{{pragma}}"
{% endfor -%}
{% endif -%}
#elif defined(_MSC_VER)
{% for name in defines -%}
#define {{name}} {{defines_table["msvc"][name]}}
{% endfor -%}
#else
{% for name in defines -%}
#define {{name}} {{defines_table["other"][name]}}
{% endfor -%}
#endif

/* section */

#if defined(__NVCC__)
{% if pragmas["nvcc"] -%}
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic pop
#else
#pragma pop
#endif // defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
{% endif -%}
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
{% if pragmas["clang_dev"] -%}
#pragma clang diagnostic pop
{% endif -%}
#else
{% if pragmas["clang_host"] -%}
// clang compiling CUDA code in host mode
#pragma clang diagnostic pop
{% endif -%}
#endif // defined(__CUDA_ARCH__)
#endif // defined(__CUDA__)
#elif defined(__GNUC__) or defined(__GNUG__)
{% if pragmas["gcc"] -%}
#pragma GCC diagnostic pop
{% endif -%}
#endif
