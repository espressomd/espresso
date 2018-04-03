import espressomd
import espressomd.checkpointing
import espressomd.virtual_sites

checkpoint = espressomd.checkpointing.Checkpointing(checkpoint_id="mycheckpoint", checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
checkpoint.register("system")
system.cell_system.skin = 0.4
system.time_step = 0.01
system.min_global_cut = 2.0
checkpoint.register("system.time_step")
checkpoint.register("system.cell_system")
checkpoint.register("system.min_global_cut")

system.part.add(pos=[1.0]*3)
system.part.add(pos=[1.0, 1.0, 2.0])

if espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']):
    system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative(have_velocity = True,
                                            have_quaternion = True)
    system.part[1].vs_auto_relate_to(0)
    checkpoint.register("system.virtual_sites")

checkpoint.register("system.part")
checkpoint.save(0)
