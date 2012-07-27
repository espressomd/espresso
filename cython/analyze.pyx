# For C-extern Analysis
cimport c_analyze

def mindist(p1 = 0, p2 = 0):
# p1,p2 should be checked for being integer list here...
	
	if p1 == 0 and p2 == 0:
		return c_analyze.mindist(NULL,NULL)
	else:
		print 'error: mindist is not correctly implemented yet\n'
#		replace "0" by "PY_ERR" or sth equal.
		return 0


def energy():
	value = 0.0
