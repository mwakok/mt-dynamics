# Microtubule dynamics simulation
Simulation of microtubule growth dynamics in the presence and absence of physical obstacles.

Work done by:
- Simulation code: Florian Huber and Maurits Kok
- Experiments run by: Maurits Kok and Svenja-Marei Kalisch
- Marileen Dogterom

## Simulation code
[![Liscence](https://img.shields.io/github/license/florian-huber/mtdynamics)](https://github.com/florian-huber/mtdynamics)
[![PyPI](https://img.shields.io/pypi/v/mtdynamics)](https://pypi.org/project/mtdynamics/0.1.0/)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)

The full simulation code used to produce the results as published in [coming soon] is provided.\ in this repository.  
The code consits of the main simulation code (**simulation_main.py** and **simulation_functions.py**). Simulation parameters are specified in **simulation_parameters.py**. Functions for plotting the results are provided in **plotting_functions.py**.

### Requirements
Python version 3.6 or higher.

### Installation
If you work with Anaconda you choose to create an own environment for ``mtdynamics`` bu running the following commands:

```
# install mtdynamics in a new virtual environment to avoid dependency clashes
conda create --name mtdynamics python=3.8
conda activate mtdynamics
pip install mtdynamics
```

Or simply install mtdynamics in your already existing environment:
```
pip install mtdynamics
```

## Runing the simulation
Jupyter notebook(s) are provided to illustrate how to run the simlation. They can be found in the folder ``\notebooks``.
