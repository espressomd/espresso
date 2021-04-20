#!/bin/bash

apt install -y $1
python3 -c 'import espressomd; system=espressomd.System(box_l=[1,1,1])'

