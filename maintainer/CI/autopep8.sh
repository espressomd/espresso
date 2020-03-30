#!/bin/bash

AUTOPEP8_VER=1.3.4
PYCODESTYLE_VER=2.3.1

if hash autopep8 2>/dev/null; then
    AUTOPEP8="$(which autopep8)"
else
    echo "No autopep8 found."
    exit 2
fi

if ! "${AUTOPEP8}" --version 2>&1 | grep -qFo "autopep8 ${AUTOPEP8_VER} (pycodestyle: ${PYCODESTYLE_VER})"; then
    echo "Could not find autopep8 ${AUTOPEP8_VER} with pycodestyle ${PYCODESTYLE_VER}"
    echo "${AUTOPEP8} is $(${AUTOPEP8} --version 2>&1)"
    exit 2
fi

${AUTOPEP8} "$@"
