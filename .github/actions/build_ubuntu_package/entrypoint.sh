#!/bin/bash
set -e -x -v
apt build-dep -y . &&
dpkg-buildpackage -uc -us &&
cp ../*.deb /github/workspace &&
package_file=`ls /github/workspace/*.deb` 
echo "::set-output name=package_file::${package_file}" 

