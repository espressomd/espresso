def AddBending(system, kb ):

    # currently only works for ONE SINGLE soft object

    # angles
    from espressomd.interactions import IBM_Tribend 
    with open("tables/softAngles", "r") as fp:
        numAngles = int( fp.readline() )
        print "Found " + str(numAngles) + " angles"    
        # actual add
        for i in range(0, numAngles):
            line = str.split( fp.readline() )
            id1 = int( line[0] )
            id2 = int( line[1] )
            id3 = int( line[2] )
            id4 = int( line[3] )            
            tribend = IBM_Tribend(ind1= id1, ind2=id2, ind3=id3, ind4=id4, kb=kb, refShape = "initial")
            system.bonded_inter.add(tribend)
            system.part[id1].add_bond((tribend, id2, id3, id4))
