#!/bin/sh
# Copyright (C) 2011 Olaf Lenz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Shell script to find out the version of ESPResSo 
#
SRCDIR=$1

test -z $SRDIR && SRCDIR=.
if ! test -d $SRCDIR; then
    echo -n "unknown"
    echo "ERROR: $SRCDIR is not a directory!" > /dev/stderr
    exit 1
fi
if ! cd $SRCDIR; then
    echo -n "unknown"
    echo "ERROR: Can't cd to $SRCDIR!" > /dev/stderr
    exit 1
fi

BASEVERSIONFILE=$SRCDIR/baseversion.txt

# try to use git describe --dirty
if ! VERSION=`git describe --dirty --match=?\.? 2> /dev/null`; then
    # try to use git without --dirty
    if ! VERSION=`git describe --match=?\.? 2> /dev/null`; then
	# otherwise use the baseversionfile
	if ! test -f $BASEVERSIONFILE; then
	    echo -n "unknown"
	    echo "ERROR: Can't find $BASEVERSIONFILE!" > /dev/stderr
	    exit 1
	else
	    VERSION=`cat $BASEVERSIONFILE`
	fi
    fi
fi
echo $VERSION | tr -d '\n'
