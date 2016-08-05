#!/bin/sh

cd `dirname $0`

libtoolize
autoreconf -iv -Wall
autoreconf -fv -Wall
