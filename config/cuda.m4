dnl CUDA support for libtool

AC_DEFUN([AC_PROG_CUDA],
      [LT_LANG(CUDA)])

AC_DEFUN([LT_PROG_CUDA],
    [AC_ARG_VAR(NVCC,[NVIDIA CUDA compiler command])
     AC_ARG_VAR(NVCCFLAGS,[special compiler flags for the NVIDIA CUDA compiler])

     AC_PATH_PROG(NVCC, nvcc, no, [$PATH:$cuda_path/bin])

     # MAC nvcc stays 32 bit, even if the rest is 64 bit
     case $target in
         x86_64-apple-darwin*)
                 NVCCFLAGS="$NVCCFLAGS -m64";;
     esac
])

# _LT_LANG_CUDA_CONFIG([TAG])
# --------------------------
# Analogue to _LT_LANG_GCJ_CONFIG for CUDA
AC_DEFUN([_LT_LANG_CUDA_CONFIG],
  [AC_REQUIRE([LT_PROG_CUDA])
   AC_LANG_PUSH(C++)

   # CUDA file extensions
   ac_ext=cu
   objext=o
   _LT_TAGVAR(objext, $1)=$objext

   # Code to be used in simple compile tests
   lt_simple_compile_test_code="static __device__ __constant__ int var;"

   # Code to be used in simple link tests
   lt_simple_link_test_code="#include <cuda.h>
                             int main() { cudaGetDevice(0); }"

   # ltmain only uses $CC for tagged configurations so make sure $CC is set.
   _LT_TAG_COMPILER

   # save warnings/boilerplate of simple test code
   _LT_COMPILER_BOILERPLATE
   _LT_LINKER_BOILERPLATE

   # Allow CC to be a program name with arguments.
   lt_save_CC="$CC"
   lt_save_GCC=$GCC

   # nvcc interface is not gcc-like (but can steer gcc)
   GCC=no
   CC=$NVCC
   compiler=$CC
   _LT_TAGVAR(compiler, $1)=$CC
   _LT_TAGVAR(LD, $1)="$LD"
   _LT_CC_BASENAME([$compiler])

   # CUDA did not exist at the time GCC didn't implicitly link libc in.
   _LT_TAGVAR(archive_cmds_need_lc, $1)=no
   _LT_TAGVAR(old_archive_cmds, $1)=$old_archive_cmds
   _LT_TAGVAR(reload_flag, $1)=$reload_flag
   _LT_TAGVAR(reload_cmds, $1)=$reload_cmds

   ## CAVEAT EMPTOR:
   ## There is no encapsulation within the following macros, do not change
   ## the running order or otherwise move them around unless you know exactly
   ## what you are doing...
   if test -n "$compiler"; then
       _LT_COMPILER_NO_RTTI($1)
       # building shared with nvcc not there in libtool
       _LT_TAGVAR(lt_prog_compiler_wl, $1)='-Xlinker '
       _LT_TAGVAR(lt_prog_compiler_static, $1)='-Xcompiler -static'
       _LT_TAGVAR(lt_prog_compiler_pic, $1)='-Xcompiler -fPIC'
       _LT_COMPILER_C_O($1)
       _LT_COMPILER_FILE_LOCKS($1)
       _LT_LINKER_SHLIBS($1)
       _LT_LINKER_HARDCODE_LIBPATH($1)

       _LT_CONFIG($1)
   fi

   AC_LANG_POP

   GCC=$lt_save_GCC
   CC="$lt_save_CC"
])
