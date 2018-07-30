#!/bin/bash

for T in *.py;
do
	if grep -q "$T" CMakeLists.txt; then
		continue;
	else
		GIT_STATUS=$(git status --porcelain -- "$T")
		if [[ $GIT_STATUS == ??* ]]; then
			echo "File '$T' is not tracked."
			continue;
		else
			echo "File '$T' is missing in CMakeLists.txt.";			
		fi
	fi
done


