add_library(Bonds SHARED PairBond.cpp Fene.cpp Harmonic.cpp
HarmonicDumbbell.cpp BondedCoulomb.cpp Quartic.cpp SubtLj.cpp Umbrella.cpp
ThreeParticleBond.cpp IbmTriel.cpp AngleHarmonic.cpp AngleCosine.cpp
AngleCosSquare.cpp AngleDist.cpp FourParticleBond.cpp HydrogenBond.cpp
MembraneCollision.cpp OifLocalForces.cpp Dihedral.cpp
TabulatedBondLength.cpp TabulatedBondAngle.cpp TabulatedBondDihedral.cpp
OverlapBondLength.cpp OverlapBondAngle.cpp OverlapBondDihedral.cpp
EightParticleBond.cpp TwistStack.cpp IbmTribend.cpp VirtualBond.cpp
BondContainer.cpp ThreeParticlePressureBond.cpp OifGlobalForces.cpp IbmVolumeConservation.cpp)
install(TARGETS Bonds LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Bonds PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Bonds EspressoConfig)
set_target_properties(Bonds PROPERTIES MACOSX_RPATH TRUE)
