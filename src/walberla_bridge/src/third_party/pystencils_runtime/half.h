/*
Copyright 2023, Markus Holzer.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


/// Half precision support. Experimental. Use carefully.
///
/// This feature is experimental, since it strictly depends on the underlying architecture and compiler support.
/// On x86 architectures, what you can expect is that the data format is supported natively only for storage and
/// interchange. Arithmetic operations will likely involve casting to fp32 (C++ float) and truncation to fp16.
/// Only bandwidth bound code may therefore benefit. None of this is guaranteed, and may change in the future.

/// Clang version must be 15 or higher for x86 half precision support.
/// GCC version must be 12 or higher for x86 half precision support.
/// Also support seems to require SSE, so ensure that respective instruction sets are enabled.
/// See
///   https://clang.llvm.org/docs/LanguageExtensions.html#half-precision-floating-point
///   https://gcc.gnu.org/onlinedocs/gcc/Half-Precision.html
/// for more information.

#pragma once
using half    = _Float16;
