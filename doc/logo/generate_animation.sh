#!/bin/sh

if ! test -d cup-animation; then
    echo "Directory cup-animation doesn't exist!"
    exit 1
fi

echo "GENERATING ANIMATED GIF..."

echo "  Combining logo box with cup in all frames..."
for i in `seq 24 49`; do
    convert -page 500x500 -gravity NorthWest logo-template.png cup-animation/$i.png -geometry +113+0 -composite cup-animation/logo-$i.png
done

echo "  Generating animation..."
for i in `seq 25 49; seq 48 -1 25`; do
    sequence="$sequence cup-animation/logo-$i.png"
done
convert \
    -dispose None \
    -delay 0 \
    cup-animation/logo-24.png \
    -dispose previous \
    -delay 5 \
    $sequence \
    -loop 0 \
    -coalesce \
    cup-animation/logo-animated-500.miff

echo "  Generating different sizes..."
convert -layers Optimize cup-animation/logo-animated-500.miff logo-animated-500.gif
convert -size 500x500 -layers Optimize cup-animation/logo-animated-500.miff -resize 200x200 logo-animated-200.gif
convert -size 500x500 -layers Optimize cup-animation/logo-animated-500.miff -resize 100x100 logo-animated-100.gif
ls -l logo-animated-500.gif logo-animated-200.gif logo-animated-100.gif

echo "Finished."
