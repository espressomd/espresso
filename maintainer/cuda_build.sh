#!/usr/bin/env bash

nvidia-smi
srcdir="${CI_PROJECT_DIR}" maintainer/CI/build_cmake.sh

