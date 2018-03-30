from espressomd import script_interface

oif_l_force = script_interface.PScriptInterface(name='Bond::OifLocalForces', policy='GLOBAL',
                                                phi0=1., kb=1., r0=1., ks=1.,kslin=1.,A01=1.,
                                                A02=1.,kal=1.)

membrane = script_interface.PScriptInterface(name='Bond::MembraneCollision', policy='GLOBAL')

dihedral_bond = script_interface.PScriptInterface(name='Bond::Dihedral',
                                                  policy='GLOBAL', mult=1.,bend=1.,phase=1.)

hydrogen_bond = script_interface.PScriptInterface(name='Bond::HydrogenBond',
                                                  policy='GLOBAL', r0=1.,alpha=1., E0=1., kd=1.,
                                                  sigma1=1., sigma2=1., psi10=1., psi20=1.,
                                                  E0sb=1., r0sb=1., alphasb=1., f2=1., f3=1.)

ibm_triel = script_interface.PScriptInterface(name='Bond::IbmTriel',
                                              policy='GLOBAL', a1=1.,a2=1.,b1=1,b2=1.,l0=1.,
                                              lp0=1., sinPhi0=1., cosPhi0=1., area0=2.,maxdist=2.,
                                              elasticLaw=3.,k1=2.,k2=1.)

angle_dist = script_interface.PScriptInterface(name='Bond::AngleDist',
                                               policy='GLOBAL', bend=1.,phimin=2., distmin=2.,
                                               phimax=3., distmax=3.)

angle_cossquare = script_interface.PScriptInterface(name='Bond::AngleCosSquare',
                                                     policy='GLOBAL', bend=1.,phi0=2.)

angle_cosine = script_interface.PScriptInterface(name='Bond::AngleCosine',
                                                     policy='GLOBAL', bend=1.,phi0=2.)

angle_harmonic = script_interface.PScriptInterface(name='Bond::AngleHarmonic',
                                                     policy='GLOBAL', bend=1.,phi0=2.)

tab_pot_dihedral = script_interface.PScriptInterface(name='Bond::TabulatedBondDihedral',
                                                     policy='GLOBAL', min=1., max=1.,
                                                     energy=[0.0,0.1,0.2,0.3], force=[1.,2.,3.,4.])

tab_pot_angle = script_interface.PScriptInterface(name='Bond::TabulatedBondAngle',
                                                  policy='GLOBAL', min=1., max=1.,
                                                  energy=[0.0,0.1,0.2,0.3], force=[1.,2.,3.,4.])

overlap_dihedral = script_interface.PScriptInterface(name='Bond::OverlapBondDihedral',
                                                     policy='GLOBAL', filename="test", maxval=1.,
                                                     noverlaps=2, para_a=[1.,2.],para_b=[1.,2.],
                                                     para_c=[1.,2.])

overlap_angle = script_interface.PScriptInterface(name='Bond::OverlapBondAngle',
                                                  policy='GLOBAL', filename="test", maxval=1.,
                                                  noverlaps=2, para_a=[1.,2.],para_b=[1.,2.],
                                                  para_c=[1.,2.])

overlap_length = script_interface.PScriptInterface(name='Bond::OverlapBondLength',
                                                   policy='GLOBAL', filename="test", maxval=1.,
                                                   noverlaps=2, para_a=[1.,2.],para_b=[1.,2.],
                                                   para_c=[1.,2.])

tab_pot_length = script_interface.PScriptInterface(name='Bond::TabulatedBondLength',
                                                   policy='GLOBAL', min=1., max=1.,
                                                   energy=[0.0,0.1,0.2,0.3], force=[1.,2.,3.,4.])

umbrella = script_interface.PScriptInterface(name='Bond::Umbrella', policy='GLOBAL', k=1.,
                                             dir=1,r=2.)

subtlj = script_interface.PScriptInterface(name='Bond::SubtLj', policy='GLOBAL')

quartic = script_interface.PScriptInterface(name='Bond::Quartic', policy='GLOBAL',
                                            k0=1., k1=2., r=1.,r_cut=10.)

bonded_coulomb_p3msr = script_interface.PScriptInterface(name='Bond::BondedCoulombP3MSR',
                                                         policy='GLOBAL', q1q2=2.0)

bonded_coulomb = script_interface.PScriptInterface(name='Bond::BondedCoulomb',
                                                   policy='GLOBAL', prefactor=1.)

harmonic = script_interface.PScriptInterface(name='Bond::Harmonic',
                                             policy='GLOBAL', k=1., r_0=1., r_cut=1.)

harmonic_dumbbell = script_interface.PScriptInterface(name='Bond::HarmonicDumbbell',
                                                      policy='GLOBAL', k1=1., k2=1., r_0=1.,
                                                      r_cut=1.)

fene = script_interface.PScriptInterface(name='Bond::Fene', policy='GLOBAL', r0=1., dr_max=2., k=1.)

bonds = script_interface.PScriptInterface(name='Bond::Bonds', policy='GLOBAL')

bonds.call_method('add', object=oif_l_force)
bonds.call_method('remove', object=oif_l_force)

bonds.call_method('add', object=membrane)
bonds.call_method('remove', object=membrane)

bonds.call_method('add', object=dihedral_bond)
bonds.call_method('remove', object=dihedral)

bonds.call_method('add', object=hydrogen_bond)
bonds.call_method('remove', object=hydrogen_bond)

bonds.call_method('add', object=ibm_triel)
bonds.call_method('remove', object=ibm_triel)

bonds.call_method('add', object=angle_dist)
bonds.call_method('remove', object=angle_dist)

bonds.call_method('add', object=angle_cossquare)
bonds.call_method('remove', object=angle_cossquare)

bonds.call_method('add', object=angle_cosine)
bonds.call_method('remove', object=angle_cosine)

bonds.call_method('add', object=angle_harmonic)
bonds.call_method('remove', object=angle_harmonic)


bonds.call_method('add', object=tab_pot_dihedral)
bonds.call_method('remove', object=tab_pot_dihedral)

bonds.call_method('add', object=tab_pot_angle)
bonds.call_method('remove', object=tab_pot_angle)

bonds.call_method('add', object=overlap_dihedral)
bonds.call_method('remove', object=overlap_dihedral)

bonds.call_method('add', object=overlap_angle)
bonds.call_method('remove', object=overlap_angle)

bonds.call_method('add', object=overlap_length)
bonds.call_method('remove', object=overlap_length)

bonds.call_method('add', object=tab_pot_length)
bonds.call_method('remove', object=tab_pot_length)

bonds.call_method('add', object=umbrella)
bonds.call_method('remove', object=umbrella)

bonds.call_method('add', object=subtlj)
bonds.call_method('remove', object=subtlj)

bonds.call_method('add', object=quartic)
bonds.call_method('remove', object=quartic)

bonds.call_method('add', object=bonded_coulomb_p3msr)
bonds.call_method('remove', object=bonded_coulomb_p3msr)

bonds.call_method('add', object=bonded_coulomb)
bonds.call_method('remove', object=bonded_coulomb)

bonds.call_method('add', object=harmonic_dumbbell)
bonds.call_method('remove', object=harmonic_dumbbell)

bonds.call_method('add', object=harmonic)
bonds.call_method('remove', object=harmonic)

bonds.call_method('add', object=fene)
bonds.call_method('remove', object=fene)
