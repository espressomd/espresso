
def changeVolume(dNew, dir="xyz"):
    if dNew<0:
        raise ValueError("No negative lengths")
    if dir=="xyz":
        dNew=dNew**(1./3.)
        rescale_boxl(3, dNew)
    elif dir=="x":
        rescale_boxl(0, dNew)
    elif dir=="y":
        rescale_boxl(1, dNew)
    elif dir=="z":
        rescale_boxl(2, dNew)
    else:
        raise ValueError('Usage: changeVolume { <V_new> | <L_new> { "x" | "y" | "z" | "xyz" } }')


    
