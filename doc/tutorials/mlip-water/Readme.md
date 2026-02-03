# Integrating MLIPs with ESPResSo: A Tutorial on Simulating Water

## Environment setup

To install Python packages and download the training data:

```sh
# create virtual environment
python -m venv venv
source venv/bin/activate
# install dependencies
pip install .
# patch dependencies
sed -i 's/periodic=np.any(atoms.pbc)/periodic=bool(np.any(atoms.pbc))/' $VIRTUAL_ENV/lib/python3.12/site-packages/apax/md/ase_calc.py
# build and install Packmol
export PACKMOL_VER=20.15.1
wget https://github.com/m3g/packmol/archive/refs/tags/v${PACKMOL_VER}.tar.gz
tar -xzf v${PACKMOL_VER}.tar.gz
rm v${PACKMOL_VER}.tar.gz
cd packmol-${PACKMOL_VER}
make -j 8
bash ./configure
make -j 8
cp ./packmol "${VIRTUAL_ENV}/bin/"
cd ..
rm -r packmol-${PACKMOL_VER}
# download training data
dvc pull
# launch JupyterLab inside the environment
PYTHONPATH=$VIRTUAL_ENV/lib/python3.12/site-packages ipypresso lab
```

## Run the tutorial

The tutorial is split into 3 parts:

### Part 1: TIP4P Water

Here, you'll learn how to run a classical MD simulation of water using the TIP4P potential.

### Part 2: MLIP Models

Here, you'll learn the fundamentals of using MLIP with ASE. You'll load and evaluate pre-trained MLIPs model for water.

### Part 3: MLIP with ESPResSo

Here, you'll learn how to set up and run an MD simulation of water using MLIP with ESPResSo.
