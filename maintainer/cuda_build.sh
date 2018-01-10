#!/usr/bin/env bash

nvidia-smi
srcdir="${CI_PROJECT_DIR}" build_procs=4 maintainer/travis/build_cmake.sh

