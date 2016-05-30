from subprocess import call
from sys import argv

skip_code = 42

returnvalue = call(argv[1:])

if returnvalue == skip_code:
    print "Skipped test '%s'" % " ".join(argv[1:])
    exit(0)
else:
    exit(returnvalue)

