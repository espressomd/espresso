def AddVolCons(system, kV ):

    # currently only works for ONE SINGLE soft object
    
    # make interaction
    from espressomd.interactions import IBM_VolCons
    volCons = IBM_VolCons(softID=1, kappaV=kV)
    system.bonded_inter.add(volCons)
    
    # loop over particles and add
    for i in range(len(system.part)):
        system.part[i].add_bond((volCons,))
