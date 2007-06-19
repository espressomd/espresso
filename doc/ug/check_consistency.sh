#!/bin/bash

#################### check documentation of configuration switches

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

doc_switches=`awk 'BEGIN { tag="feature{"; l=length(tag) }
    /feature/ {
	s=$0;
	while(1) {
	    pos=index(s, "feature{");
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
