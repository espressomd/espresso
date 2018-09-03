import unittest as ut

from ek_eof_one_species_base import ek_eof_one_species
from ek_eof_one_species_base import params_base

params_z = dict([
    ('box_x', 3.0),
    ('box_y', params_base['width'] + 2 * params_base['padding']),
    ('box_z', 3.0),
    ('ext_force_density', [0.0, 0.0, params_base['force']]),
    ('wall_normal_1', [0, 1, 0]),
    ('wall_normal_2', [0, -1, 0]),
    ('periodic_dirs', (0, 2)),
    ('non_periodic_dir', 1),
    ('n_roll_index', 2),
    ('calculated_pressure_xy', 0.0),
    ('calculated_pressure_xz', 0.0)
])


class eof_z(ek_eof_one_species):
    def test(self):
        self.run_test(params_z)


if __name__ == "__main__":
    ut.main()
