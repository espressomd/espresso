# Check whether all features used in the code are defined
import featuredefs

if len(sys.argv) != 2:
    print "Usage: %s FILE" % sys.argv[0]
    exit(2)

defs = featuredefs.defs(sys.argv[1])

