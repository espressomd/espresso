/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

int ObservableAverage::actual_update() {
  observable_average_container *data =
      (observable_average_container *)container;
  data->n_sweeps++;
  int error = data->reference_observable->calculate();
  if (error != 0)
    return 1;
  double factor = 1 / (double)data->n_sweeps;
  for (int i = 0; i < n; i++) {
    last_value[i] = (1 - factor) * last_value[i] +
                    factor * data->reference_observable->last_value[i];
  }
  return 0;
}

int ObservableAverage::reset() {
  observable_average_container *data =
      (observable_average_container *)container;
  data->n_sweeps = 0;
  int error = data->reference_observable->calculate();
  for (int i = 0; i < n; i++) {
    last_value[i] = 0;
  }
  return error;
}
