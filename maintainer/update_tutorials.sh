#!/bin/bash

# Build tutorials and copy and add
# the new pdf's into place.
# To be run from the build directory.

make tutorials

for F in $(find -type f -name "*.pdf" | xargs); do 
	TARGET="../$F";
        echo "Updating ${TARGET}"
	cp $F ${TARGET}; 
	git add $TARGET; 
done

