from espressomd.interactions import BondedCoulombP3MSRBond

#Dict with drude type infos
drude_dict={}
#Lists with unique drude and core types
core_type_list=[]
drude_type_list=[]
#Get core id from drude id
core_id_from_drude_id={}

def add_drude_particle_to_core(system, p_core, drude_bond, id_drude, type_drude, alpha, mass_drude, coulomb_prefactor, thole_damping = 2.6):
    
    k = drude_bond.params["k"]
    
    k=drude_bond.params["k"]
    q_drude =  -1.0 * pow(k * alpha / coulomb_prefactor, 0.5)

    system.part.add(id=id_drude, pos=p_core.pos, type = type_drude,  q = q_drude, mass = mass_drude, temp = 0, gamma = 0)
    #print("Adding to core", p_core.id, "drude id", id_drude, "  pol", alpha, "  core charge", p_core.q, "->", p_core.q-q_drude, "   drude charge", q_drude)

    p_core.q -= q_drude
    p_core.mass -= mass_drude   
    p_core.add_bond((drude_bond, id_drude))
    p_core.temp = 0
    p_core.gamma = 0

    #THOLE
    #See LAMMPS paper on Drude (DOI 10.1021/acs.jcim.5b00612)
    
    #Drude particles with different drude charges(or alphas)/thole_damping have to have different types for thole:
    #global drude_dict
    if type_drude in drude_dict and not (drude_dict[type_drude]["q"] == q_drude and drude_dict[type_drude]["thole_damping"] == thole_damping):
        print("ERROR: Drude particles with different drude charges have to have different types for thole..Aborting")
        exit()
        
    core_id_from_drude_id[id_drude] = p_core.id


    #Add new thole nonbonded interaction for D-D, D-C, C-C for all existing drude types if this type is seen for the first time
    if not type_drude in drude_dict:
    
        #Bookkepping of q, alphas and damping parameter
        drude_dict[type_drude] = {}
        drude_dict[type_drude]["q"] = q_drude
        drude_dict[type_drude]["alpha"] = alpha
        drude_dict[type_drude]["thole_damping"] = thole_damping
        drude_dict[type_drude]["core_type"] = p_core.type
        #Save same information to get access to the parameters via core types
        drude_dict[p_core.type] = {}
        drude_dict[p_core.type]["q"] = -q_drude
        drude_dict[p_core.type]["alpha"] = alpha
        drude_dict[p_core.type]["thole_damping"] = thole_damping
        drude_dict[p_core.type]["drude_type"] = type_drude

    #Collect unique drude types
    if not type_drude in drude_type_list:
        drude_type_list.append(type_drude)
    
    #Collect unique core types
    if not p_core.type in core_type_list:
        core_type_list.append(p_core.type)


def add_thole_pair_damping(system, t1,t2):
    qq = drude_dict[t1]["q"] * drude_dict[t2]["q"]
    s = 0.5 * (drude_dict[t1]["thole_damping"] + drude_dict[t2]["thole_damping"]) / (drude_dict[t1]["alpha"] * drude_dict[t2]["alpha"])**(1.0/6.0) 
    system.non_bonded_inter[t1,t2].thole.set_params(scaling_coeff=s, q1q2 = qq)
    #print("Added THOLE for types", t1,"<->", t2, "S",s, "q1q2",qq)

def add_all_thole(system):
    #drude <-> drude
    for i in range(len(drude_type_list)):
        for j in range(i,len(drude_type_list)):
            add_thole_pair_damping(system, drude_type_list[i], drude_type_list[j])
    #core <-> core
    for i in range(len(core_type_list)):
        for j in range(i,len(core_type_list)):
            add_thole_pair_damping(system, core_type_list[i], core_type_list[j])
    #drude <-> core
    for i in drude_type_list:
        for j in core_type_list:
            add_thole_pair_damping(system, i, j)

def setup_intramol_exclusion_bonds(system, mol_drude_types, mol_core_types, mol_core_partial_charges):

        #All drude types need...
        for td in mol_drude_types:
            drude_dict[td]["subtr_p3m_sr_bonds"]={}

            #...p3m sr exclusion bond with other partial core charges...
            for tc, qp in zip(mol_core_types, mol_core_partial_charges):
                #...excluding the drude core partner
                if drude_dict[td]["core_type"] != tc:
                    qd = drude_dict[td]["q"] #Drude charge
                    subtr_p3m_sr_bond = BondedCoulombP3MSRBond(q1q2 = -qd*qp)
                    system.bonded_inter.add(subtr_p3m_sr_bond)
                    drude_dict[td]["subtr_p3m_sr_bonds"][tc]=subtr_p3m_sr_bond
                    #print("Added intramolecular exclusion", subtr_p3m_sr_bond, "for drude",  qd, "<-> core", qp, "to system") 
        
def add_intramol_exclusion_bonds(system, drude_ids, core_ids):

    for drude_id in drude_ids:
        for core_id in core_ids:
            if core_id_from_drude_id[drude_id] != core_id:
                pd = system.part[drude_id]
                pc = system.part[core_id]
                bond = drude_dict[pd.type]["subtr_p3m_sr_bonds"][pc.type]
                pd.add_bond((bond, core_id))
                #print("Added subtr_p3m_sr bond", bond, "between ids", drude_id, "and", core_id)

