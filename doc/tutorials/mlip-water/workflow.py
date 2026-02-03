#
# Copyright (C) 2025-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import yaml
import copy
import pathlib
import ipsuite as ips
from apax.nodes import Apax, ApaxBatchPrediction


class Workflow:

    def __init__(self):
        self.reference = yaml.safe_load(
            pathlib.Path("config/apax.yaml").read_text())
        self.project = ips.Project()
        self.configured = False

    def configure(self):
        with self.project:
            train = ips.AddDataH5MD(file="data/cosmo_water_train.h5")
            val = ips.AddDataH5MD(file="data/cosmo_water_val.h5")
            test = ips.AddDataH5MD(file="data/cosmo_water_test.h5")

        for r_max in [2, 3, 4, 5, 5.5, 6]:
            with self.project.group("r_max", f"{r_max}".replace(".", "_")):
                config = copy.deepcopy(self.reference)
                config["model"]["basis"]["r_max"] = r_max
                pth = pathlib.Path(f"config/apax-r_max-{r_max}.yaml")
                pth.write_text(yaml.dump(config))
                model = Apax(
                    config=pth.as_posix(),
                    data=train.frames,
                    validation_data=val.frames,
                )
                test_eval = ApaxBatchPrediction(data=test.frames, model=model)
                ips.PredictionMetrics(x=test_eval.frames, y=test.frames)

        for nn in [(16, 16), (32, 32), (64, 64), (128, 128)]:
            with self.project.group("nn", f"{nn[0]}-{nn[1]}"):
                config = copy.deepcopy(self.reference)
                config["model"]["nn"] = list(nn)
                pth = pathlib.Path(f"config/apax-nn-{nn[0]}-{nn[1]}.yaml")
                pth.write_text(yaml.dump(config))
                model = Apax(
                    config=pth.as_posix(),
                    data=train.frames,
                    validation_data=val.frames,
                )
                test_eval = ApaxBatchPrediction(data=test.frames, model=model)
                ips.PredictionMetrics(x=test_eval.frames, y=test.frames)

        for n_basis in [4, 8, 16]:
            with self.project.group("n_basis", f"{n_basis}"):
                config = copy.deepcopy(self.reference)
                config["model"]["basis"]["n_basis"] = n_basis
                pth = pathlib.Path(f"config/apax-n_basis-{n_basis}.yaml")
                pth.write_text(yaml.dump(config))
                model = Apax(
                    config=pth.as_posix(),
                    data=train.frames,
                    validation_data=val.frames,
                )
                test_eval = ApaxBatchPrediction(data=test.frames, model=model)
                ips.PredictionMetrics(x=test_eval.frames, y=test.frames)

        for n_radial in [5, 6, 7]:
            with self.project.group("n_radial", f"{n_radial}"):
                config = copy.deepcopy(self.reference)
                config["model"]["n_radial"] = n_radial
                pth = pathlib.Path(f"config/apax-n_radial-{n_radial}.yaml")
                pth.write_text(yaml.dump(config))
                model = Apax(
                    config=pth.as_posix(),
                    data=train.frames,
                    validation_data=val.frames,
                )
                test_eval = ApaxBatchPrediction(data=test.frames, model=model)
                ips.PredictionMetrics(x=test_eval.frames, y=test.frames)

        self.configured = True

    def train(self):
        assert self.configured, "project hasn't been configured yet"
        self.project.build()
