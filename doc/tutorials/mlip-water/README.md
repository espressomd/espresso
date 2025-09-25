# Integration MLIPs with ESPResSo: A Tutorial on Simulating Water

## Setup (via UV)

1. Install [UV](https://astral.sh/uv/) if you haven't already `curl -LsSf https://astral.sh/uv/install.sh | sh`
2. Create the virtual environment and download the data and model files

```bash
uv sync --frozen --no-cache # --no-cache is optional to avoid caching into your home directory
source .venv/bin/activate # activate venv
dvc pull # download the model and data files
```

## Setup (via pip)
You can also install the tutorial via pip.
This will not use the `uv.lock` file and dependencies version can be different.

```bash
pip install .
dvc pull
```

## Run the tutorial

The tutorial is split into 3 parts:
### TIP4P Water
In the first part, you'll learn how to run a classical MD simulation of water using the TIP4P potential.

### MLIP Models
In the second part, you'll learn the fundamentals of using MLIP with ASE. You'll load and evaluate pre-trained MLIPs model for water.

### MLIP with ESPResSo
In the third part, you'll learn how to set up and run an MD simulation of water using MLIP with ESPResSo.
