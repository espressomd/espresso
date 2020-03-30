#!/bin/bash


CLANG_FORMAT_VER=6.0
if hash clang-format_${CLANG_FORMAT_VER} 2>/dev/null; then
    CLANGFORMAT="$(which clang-format-${CLANG_FORMAT_VER})"
elif hash clang-format 2>/dev/null; then
    CLANGFORMAT="$(which clang-format)"
else
    echo "No clang-format found."
    exit 2
fi

if ! "${CLANGFORMAT}" --version | grep -qEo "version ${CLANG_FORMAT_VER}\.[0-9]+"; then
    echo "Could not find clang-format ${CLANG_FORMAT_VER}. ${CLANGFORMAT} is $(${CLANGFORMAT} --version | grep -Eo '[0-9\.]{5}' | head -n 1)."
    exit 2
fi

${CLANGFORMAT} "$@"
