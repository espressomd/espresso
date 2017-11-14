from __future__ import print_function
from sys import argv
import sys
from subprocess import CalledProcessError

from defines import Defines
import featuredefs

def damerau_levenshtein_distance(s1, s2):
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in range(-1,lenstr1+1):
        d[(i,-1)] = i+1
    for j in range(-1,lenstr2+1):
        d[(-1,j)] = j+1

    for i in range(lenstr1):
        for j in range(lenstr2):
            if s1[i] == s2[j]:
                cost = 0
            else:
                cost = 1
            d[(i,j)] = min(
                           d[(i-1,j)] + 1, # deletion
                           d[(i,j-1)] + 1, # insertion
                           d[(i-1,j-1)] + cost, # substitution
                          )
            if i and j and s1[i]==s2[j-1] and s1[i-1] == s2[j]:
                d[(i,j)] = min (d[(i,j)], d[i-2,j-2] + cost) # transposition

    return d[lenstr1-1,lenstr2-1]

def handle_unkown(f, all_features):
    match = None
    max_dist = max(2, len(f) // 2)
    for d in all_features:
        dist = damerau_levenshtein_distance(f, d)
        if dist < max_dist:
            min_dist = dist
            match = d

    if match:
        print("Unknown feature '{}', did you mean '{}'?".format(f, match))
    else:
        print("Unknown feature '{}'".format(f))

class ParseError(Exception):
    pass

def check_myconfig(compiler, feature_file, myconfig):
# This does not work on all compilers, so if the parsing fails
# we just bail out.
    try:
        my_features = Defines(compiler).defines(myconfig)
    except CalledProcessError as ex:
        print("""Skipped myconfig check because {} returned non-zero exit code {},
                  output: {}.""".format(ex.cmd, ex.returncode,ex.output.strip()))
        return

# Parse feature file
    defs = featuredefs.defs(feature_file)

    error_state = False
    for e in (my_features & defs.externals):
        error_state = True
        my_features.remove(e)
        print("External feature '{}' can not be defined in myconfig.".format(e))

    for u in (my_features - defs.features):
        error_state = True
        handle_unkown(u, defs.features)

    if error_state:
        raise ParseError("There were errors in '{}'".format(argv[3]))
    else:
        return

if __name__ == "__main__":
    try:
        check_myconfig(argv[1], argv[2], argv[3])
        sys.exit()
    except ParseError:
        sys.exit("There were errors in '{}'".format(argv[3]))
