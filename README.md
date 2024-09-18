# pseudo-hard-sphere-solvation-paper
Repository containing setup and analysis scripts used in the paper studying the high-energy water in cucurbit[8]uril.


NOTE: This repo is still a work in progress, there are still many things I will add in the near future.

The main functions that converts/creates the hard sphere model of a dummy particle is in `hard_sphere.py`.

The file `create_hard-sphere.py` creates an OpenMM System XML file containing the modified potentials after loading the AMBER files `*.prmtop` and `*.rst7`

TODO: Add a script that runs the simulations (with and without umbrella sampling restraints) and the alchemical version that was used to estiamte the hydration free energy.
