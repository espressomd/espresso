#!/bin/bash
# check documentation of configuration switches
#
# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2007,2008,2009,2010 Axel Arnold
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.# 
#

#### list comparison helper
listdiff() {
    for t in $1; do
	eval red=\${3%%$t*}
	if [ "$red" == "$3" ]; then
	    echo "$t $2" 1>&2
	fi
    done
    for t in $3; do
 	eval red=\${1%%$t*}
 	if [ "$red" == "$1" ]; then
 	    echo "$t $4" 1>&2
 	fi
    done
}

#### 1. get switches from features.tex

doc_switches=`awk 'BEGIN { tag="newfeature{"; l=length(tag) }
    /newfeature/ {
	s=$0;
	while(1) {
	    pos=index(s, "newfeature{");
	    if (pos==0) break;
	    s=substr(s,pos+l);
	    pos=index(s, "}");
	    if (pos==0) {
		print "ERROR: missing closing brace after feature in line "$0 > "/dev/stderr";
		pos=length(s)
	    };
	    tag=substr(s, 1, pos-1);
	    gsub("\\\\\\\\_", "_", tag);
	    print tag;
	    s=substr(s, pos+1);
	}
    }' features.tex | sort`

# echo "found the following documented switches: $doc_switches"

#### 1. get switches from myconfig-sample.h
have_switches=`awk '/\\/\\* #define/ { for(i=0;i<NF;i++) { if ($i=="#define") { print $(i+1); break} } }' ../../myconfig-sample.h | sort`

# echo "found the following possible switches: $have_switches"

#### and print difference

listdiff "$doc_switches" "documented, but not in myconfig-sample.h" "$have_switches" "not documented in features.tex"
