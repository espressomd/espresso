#!/bin/bash

apt-get install -y $1
echo $1
ls -la /github/workspace
python3 -c 'import espressomd; system=espressomd.System(box_l=[1,1,1])'
