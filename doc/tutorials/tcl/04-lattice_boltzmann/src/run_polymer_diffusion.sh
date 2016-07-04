#!/usr/bin/env sh

for i in 20 40 60 100
do
   	../espresso-3.3.0/build/Espresso polymer_diffusion.tcl $i;
done
