from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper

class AutoUpdateObservables(ScriptInterfaceHelper):
    _so_name = "Observables::AutoUpdateObservables"

    def add(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], Observable):
                observable = args[0]
            else:
                raise TypeError(
                    "Either a Observable object or key-value pairs for the parameters of a Observable object need to be passed.")
        else:
            observable = Observable(**kwargs)
        self.call_method("add", object=observable)
        return observable

    def remove(self, observable):
        self.call_method("remove", object=observable)

class Observable(ScriptInterfaceHelper):
    _so_name="Observables::Observable"
    _so_bind_methods = ("value","calculate","update","auto_write_to")



class ComForce(Observable):
    _so_name="Observables::ComForce"


class ComPosition(Observable):
    _so_name="Observables::ComPosition"


class ComVelocity(Observable):
    _so_name="Observables::ComVelocity"


class Current(Observable):
    _so_name="Observables::Current"


class DensityProfile(Observable):
    _so_name="Observables::DensityProfile"


class DipoleMoment(Observable):
    _so_name="Observables::DipoleMoment"


class FluxDensityProfile(Observable):
    _so_name="Observables::FluxDensityProfile"


class ForceDensityProfile(Observable):
    _so_name="Observables::ForceDensityProfile"


class LBVelocityProfile(Observable):
    _so_name="Observables::LBVelocityProfile"


class MagneticDipoleMoment(Observable):
    _so_name="Observables::MagneticDipoleMoment"


class ParticleAngularMomentum(Observable):
    _so_name="Observables::ParticleAngularMomentum"


class ParticleBodyAngularMomentum(Observable):
    _so_name="Observables::ParticleBodyAngularMomentum"


class ParticleBodyVelocities(Observable):
    _so_name="Observables::ParticleBodyVelocities"


class ParticleCurrent(Observable):
    _so_name="Observables::ParticleCurrent"


class ParticleForces(Observable):
    _so_name="Observables::ParticleForces"


class ParticlePositions(Observable):
    _so_name="Observables::ParticlePositions"


class ParticleVelocities(Observable):
    _so_name="Observables::ParticleVelocities"


class StressTensor(Observable):
    _so_name="Observables::StressTensor"


class StressTensorAcf(Observable):
    _so_name="Observables::StressTensorAcf"


