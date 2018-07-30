#!/usr/bin/env bash

nvidia-smi
srcdir="${CI_PROJECT_DIR}" build_procs=4 maintainer/CI/build_cmake.sh

