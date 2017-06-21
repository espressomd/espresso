from __future__ import print_function, absolute_import
from .script_interface import ScriptInterfaceHelper,script_interface_register

@script_interface_register
class AutoUpdateObservables(ScriptInterfaceHelper):
    _so_name = "Observables::AutoUpdateObservables"
    _so_creation_policy = "LOCAL"

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

@script_interface_register
class Observable(ScriptInterfaceHelper):
    _so_name="Observables::Observable"
    _so_bind_methods = ("value","calculate","update","auto_write_to")
    _so_creation_policy = "LOCAL"

@script_interface_register
class ComForce(Observable):
    _so_name="Observables::ComForce"


@script_interface_register
class ComPosition(Observable):
    _so_name="Observables::ComPosition"


@script_interface_register
class ComVelocity(Observable):
    _so_name="Observables::ComVelocity"


@script_interface_register
class Current(Observable):
    _so_name="Observables::Current"


@script_interface_register
class DensityProfile(Observable):
    _so_name="Observables::DensityProfile"


@script_interface_register
class DipoleMoment(Observable):
    _so_name="Observables::DipoleMoment"


@script_interface_register
class FluxDensityProfile(Observable):
    _so_name="Observables::FluxDensityProfile"


@script_interface_register
class ForceDensityProfile(Observable):
    _so_name="Observables::ForceDensityProfile"


@script_interface_register
class LBVelocityProfile(Observable):
    _so_name="Observables::LBVelocityProfile"


@script_interface_register
class MagneticDipoleMoment(Observable):
    _so_name="Observables::MagneticDipoleMoment"


@script_interface_register
class ParticleAngularMomentum(Observable):
    _so_name="Observables::ParticleAngularMomentum"


@script_interface_register
class ParticleBodyAngularMomentum(Observable):
    _so_name="Observables::ParticleBodyAngularMomentum"


@script_interface_register
class ParticleBodyVelocities(Observable):
    _so_name="Observables::ParticleBodyVelocities"


@script_interface_register
class ParticleCurrent(Observable):
    _so_name="Observables::ParticleCurrent"


@script_interface_register
class ParticleForces(Observable):
    _so_name="Observables::ParticleForces"


@script_interface_register
class ParticlePositions(Observable):
    _so_name="Observables::ParticlePositions"


@script_interface_register
class ParticleVelocities(Observable):
    _so_name="Observables::ParticleVelocities"


@script_interface_register
class StressTensor(Observable):
    _so_name="Observables::StressTensor"


@script_interface_register
class StressTensorAcf(Observable):
    _so_name="Observables::StressTensorAcf"


