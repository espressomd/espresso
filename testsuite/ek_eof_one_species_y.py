import unittest as ut

from ek_eof_one_species_base import ek_eof_one_species
from ek_eof_one_species_base import params_base

params_y = dict([
    ('box_x', params_base['width'] + 2 * params_base['padding']),
    ('box_y', 3.0),
    ('box_z', 3.0),
    ('ext_force_density', [0.0, params_base['force'], 0.0]),
    ('wall_normal_1', [1, 0, 0]),
    ('wall_normal_2', [-1, 0, 0]),
    ('periodic_dirs', (1, 2)),
    ('non_periodic_dir', 0),
    ('n_roll_index', 1),
    ('calculated_pressure_xz', 0.0),
    ('calculated_pressure_yz', 0.0)
])


class eof_y(ek_eof_one_species):
    def test(self):
        self.run_test(params_y)


if __name__ == "__main__":
    ut.main()
