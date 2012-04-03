#!/bin/sh

if ! test -d cup-animation; then
    echo "Directory cup-animation doesn't exist!"
    exit 1
fi

echo "GENERATING ANIMATED GIF..."
for i in `seq 24 50; seq 49 -1 24`; do
    sequence="$sequence -page 200x200+100+0 cup-animation/$i.png"
done
convert \
    -page 500x500 \
    -dispose None \
    logo-template.png \
    -set dispose previous \
    -loop 0 -delay 15 \
    $sequence \
    cup-animation-500.gif
ls -l cup-animation-500.gif

