#!/bin/bash
apt build-dep -y .
dpkg-build-package -uc -us

