## MD Simulations of pseudo-hard-sphere particle with CB8 system
There are two Python scripts in this folder:
* [simulate_z=0.py](simulate_z=0.py)
* [simulate_z=15.py](simulate_z=15.py)

The scripts creates the pseudo-hard-sphere particle with the unperturbed **CB8** model solvated in TIP3P water molecules. 

The second part of the script runs an initial energy minimization to remove any steric clashes followed by a short 200ps MD simulations. Finally, a further 5ns of MD simulations is run for the production.

Simply run the script in your local machine or cluster (I ran the MD simulations on my M2 macbook pro using OpenCL platform)

```python
python simulate_z=0.py
```
```python
python simulate_z=15.py
```

## Simulation video
Below is video of the MD simulations for $z=0$ Angstrom. 

> * The CB8 is represented as licorice
> * The pseudo-hard-sphere particle is represented as a sphere with a radius of 5 Ang.
> * The vdW spheres are the closest water molecules around the pseudo-hard-sphere particle.
> * All other water molecules are represented as lines.
> * The periodic boundary is shown in blue.

<img src="z_initial/movie.gif" width="840" height="997">

Below is video of the MD simulations for $z=15$ Angstrom.

> * The CB8 is represented as licorice
> * The pseudo-hard-sphere particle is represented as a sphere with a radius of 5 Ang.
> * The vdW spheres are the closest water molecules around the pseudo-hard-sphere particle.
> * All other water molecules are represented as lines.
> * The periodic boundary is shown in blue.