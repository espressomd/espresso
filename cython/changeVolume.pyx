
def changeVolume(dNew, dir="xyz"):
	if dir=="xyz":
		print dNew
		dNew=dNew**(1./3.)
		print dNew
		rescale_boxl(3, dNew)
	elif dir=="x":
		rescale_boxl(0, dNew)
	elif dir=="y":
		rescale_boxl(1, dNew)
	elif dir=="z":
		rescale_boxl(2, dNew)


	
