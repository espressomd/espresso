#!/bin/bash
apt build-dep -y . &&
dpkg-buildpackage -uc -us &&
ls ~/*.deb

