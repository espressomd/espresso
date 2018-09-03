import unittest as ut

from ek_eof_one_species_base import ek_eof_one_species
from ek_eof_one_species_base import params_base

params_x = dict([
    ('box_x', 3.0),
    ('box_y', 3.0),
    ('box_z', params_base['width'] + 2 * params_base['padding']),
    ('ext_force_density', [params_base['force'], 0.0, 0.0]),
    ('wall_normal_1', [0, 0, 1]),
    ('wall_normal_2', [0, 0, -1]),
    ('periodic_dirs', (0, 1)),
    ('non_periodic_dir', 2),
    ('n_roll_index', 0),
    ('calculated_pressure_xy', 0.0),
    ('calculated_pressure_yz', 0.0)
])


class eof_x(ek_eof_one_species):
    def test(self):
        self.run_test(params_x)


if __name__ == "__main__":
    ut.main()
