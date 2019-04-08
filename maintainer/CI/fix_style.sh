#!/bin/bash

cd "$(git rev-parse --show-toplevel)"

CLANG_FORMAT_VER=6.0
AUTOPEP8_VER=0.9.1

if ! git diff-index --quiet HEAD -- && [ "$1" != "-f" ]; then
    echo "Warning, your working tree is not clean. Please commit your changes."
    echo "You can also call this script with the -f flag to proceed anyway, but"
    echo "you will then be unable to revert the formatting changes later."
    exit 1
fi

CLANGFORMAT="$(which clang-format-$CLANG_FORMAT_VER)"
if [ "$CLANGFORMAT" = "" ]; then
    CLANGFORMAT="$(which clang-format)"
    if ! $CLANGFORMAT --version | grep -qEo "version $CLANG_FORMAT_VER\.[0-9]+"; then
        echo "Could not find clang-format $CLANG_FORMAT_VER. $CLANGFORMAT is $($CLANGFORMAT --version | grep -Eo '[0-9\.]{5}' | head -n 1)."
        exit 2
    fi
fi

AUTOPEP8="$(which autopep8)"
if ! $AUTOPEP8 --version | grep -qEo "autopep8 $AUTOPEP8_VER"; then
    echo "Could not find autopep8 $AUTOPEP8_VER. $AUTOPEP8 is $($AUTOPEP8 --version | grep -Eo '[0-9\.]{5}')."
    exit 2
fi

find . \( -name '*.hpp' -o -name '*.cpp' -o -name '*.cu' \) -not -path './libs/*' | xargs -r -n 5 -P 8 $CLANGFORMAT -i -style=file || exit 3
find . \( -name '*.py' -o -name '*.pyx' -o -name '*.pxd' \) -not -path './libs/*' | xargs -r -n 5 -P 8 $AUTOPEP8 --ignore=E266,W291,W293 --in-place --aggressive --aggressive || exit 3
find . -type f -executable ! -name '*.sh' ! -name '*.py' ! -name '*.sh.in' ! -name pypresso.cmakein  -not -path './.git/*' | xargs -r -n 5 -P 8 chmod -x || exit 3

if [ "$CI" != "" ]; then
    maintainer/gh_post_style_patch.py
fi

