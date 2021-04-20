#!/bin/bash
apt build-dep .
dpkg-build-package -uc -us

