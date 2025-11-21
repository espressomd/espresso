# ESPResSo in Codespaces

## Creating a new codespace

Direct link to the default codespace: https://codespaces.new/espressomd/espresso

If you need to create a custom codespace, follow these instructions:

* go to the ESPResSo repository, either the [official project](https://github.com/espressomd/espresso) or your forked project
* follow the GitHub documentation on [Creating a codespace for a repository](https://docs.github.com/en/codespaces/developing-in-a-codespace/creating-a-codespace-for-a-repository#creating-a-codespace-for-a-repository)

## Building ESPResSo

General workflow:

* load MPI module
* build ESPResSo
* create Python environment

Command lines:

```sh
# load MPI module
. /etc/profile.d/modules.sh
module load mpi
# build ESPResSo
mkdir build
cd build
cmake .. -D ESPRESSO_BUILD_WITH_FFTW=ON -D ESPRESSO_BUILD_WITH_WALBERLA=ON
make -j$(nproc)
# create Python environment
python -m venv --system-site-packages ESPResSo
. ESPResSo/bin/activate
realpath src/python > $(python -c 'import sysconfig;print(sysconfig.get_path("platlib"))')/espresso.pth
pip install -c ../requirements.txt matplotlib pint pandas tqdm
```

## Running interactive notebooks

Two options:

* run tutorials from subfolder `doc/tutorials` in the top-level directory
* run `make -j$(nproc) tutorials` to generate notebooks with hidden answers,
  and then run tutorials from subfolder `build/doc/tutorials`

The IDE needs to know which Python interpreter to use:

* open the Command Palette with F1 (or click on the gears icon)
* search for "Python: Set Project Environment"
* hit "Browse"
* paste `/workspaces/espresso/build/ESPResSo/bin/activate`

When the IDE asks for the interpreter, select the `ESPResSo` environment.
