import os
import sys

import openmm
import openmm.app as app
import openmm.unit as unit

sys.path.insert(0, "../")
from hard_sphere import make_guest_wca_repulsive

# User Options
# ----------------------------------------------------------------------- #
guest_resname = "DM1"
host_resname = "CB8"
solvent_resname = "HOH"

hard_sphere_sigma = 3.0 * unit.angstrom
hard_sphere_epsilon = 0.1 * unit.kilocalorie_per_mole
hard_sphere_radius = 5.0 * unit.angstrom

wca_repulsive = 50
wca_attractive = 49

temperature = 298.15 * unit.kelvin
pressure = 1.0 * unit.bar

work_dir = "z_final"
os.makedirs(work_dir, exist_ok=True)

# ----------------------------------------------------------------------- #
# 01 - Initial Setup
# ----------------------------------------------------------------------- #

# Load Amber Files - z=0 Ang
# ----------------------------------------------------------------------- #
prmtop = app.AmberPrmtopFile("../structures/CB8-Guest-sol.prmtop")
inpcrd = app.AmberInpcrdFile("../structures/CB8-Guest-sol-z=15.rst7")

# Get atom indices of different components in system
# ----------------------------------------------------------------------- #
guest = [
    atom.index for atom in prmtop.topology.atoms() if atom.residue.name == guest_resname
]
host = [
    atom.index for atom in prmtop.topology.atoms() if atom.residue.name == host_resname
]
solvent = [
    atom.index
    for atom in prmtop.topology.atoms()
    if atom.residue.name == solvent_resname
]

# Select Carbon atoms of CB8 to get cavity center
host_com = [
    atom.index
    for atom in prmtop.topology.atoms()
    if atom.residue.name == host_resname and atom.element.symbol == "C"
]

# Create OpenMM System
# ----------------------------------------------------------------------- #
system = prmtop.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=9.0 * unit.angstrom,
    constraints=app.HBonds,
    rigidWater=True,
)

# Set Host mass to zero (makes CB8 host rigid throughout MD simulation)
for atom in host:
    system.setParticleMass(atom, 0.0 * unit.dalton)

# Set Dummy particle to "hard-sphere"
# ----------------------------------------------------------------------- #
make_guest_wca_repulsive(
    system,
    guest,
    host,
    solvent,
    particle_radius=hard_sphere_radius,
    wca_sigma=hard_sphere_sigma,
    wca_epsilon=hard_sphere_epsilon,
    coeff_r=wca_repulsive,
    coeff_a=wca_attractive,
    force_group=10,
)

# Set Dummy particle LJ parameters to zero
# ----------------------------------------------------------------------- #
nonbonded = [
    force for force in system.getForces() if isinstance(force, openmm.NonbondedForce)
][0]
for atom in guest:
    charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
    nonbonded.setParticleParameters(atom, 0.0, sigma, 0.0)

# Apply Harmonic restraints on the Dummy Particle at z=0 Ang
# ----------------------------------------------------------------------- #
restraints = openmm.CustomCentroidBondForce(
    2,
    "(k_z/2)*(z-z0)^2 + k_xy*((x-x0)^2 + (y-y0)^2);"
    "x = x2-x1;"
    "y = y2-y1;"
    "z = z2-z1;",
)
restraints.addGlobalParameter("k_z", 25 * unit.kilocalorie_per_mole / unit.angstrom**2)
restraints.addGlobalParameter(
    "k_xy", 100 * unit.kilocalorie_per_mole / unit.angstrom**2
)
restraints.addGlobalParameter("x0", 0.0 * unit.angstrom)
restraints.addGlobalParameter("y0", 0.0 * unit.angstrom)
restraints.addGlobalParameter("z0", 0.0 * unit.angstrom)
g1 = restraints.addGroup(host_com, [1.0] * len(host_com))
g2 = restraints.addGroup(guest, [1.0] * len(guest))
restraints.addBond([g1, g2], [])
restraints.setName("Dummy Particle Restraints")
restraints.setForceGroup(10)
system.addForce(restraints)

# Save System to XML file
# ----------------------------------------------------------------------- #
with open("z_initial/system_hard-sphere.xml", "w") as file:
    file.write(openmm.XmlSerializer.serialize(system))

# ----------------------------------------------------------------------- #
# 02 - MD Simulations - 200ps Equilbration with 2.0fs timestep
# ----------------------------------------------------------------------- #
# Simulation Object
barostat = openmm.MonteCarloBarostat(pressure, temperature, 25)
system.addForce(barostat)

integrator = openmm.LangevinMiddleIntegrator(
    temperature, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
)
simulation = app.Simulation(
    prmtop.topology,
    system,
    integrator,
    openmm.Platform.getPlatformByName("OpenCL"),
)
simulation.context.setPositions(inpcrd.positions)

# Reporters
simulation.reporters.append(
    app.CheckpointReporter(f"{work_dir}/checkpoint_equilibration.xml", 25000, writeState=True)
)
simulation.reporters.append(
    app.DCDReporter(f"{work_dir}/trajectory_equilibration.dcd", 2500)
)
simulation.reporters.append(
    app.StateDataReporter(
        f"{work_dir}/reporter_equilibration.csv",
        2500,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
    )
)

# Minimization
simulation.minimizeEnergy()

# MD Simulation
simulation.step(100000)

# Write final snapshot
simulation.saveState(f"{work_dir}/checkpoint_equilibration.xml")
with open(f"{work_dir}/complex_equilibration.pdb", "w") as f:
    app.PDBFile.writeFile(
        simulation.topology,
        simulation.context.getState(getPositions=True).getPositions(),
        f,
        keepIds=True,
    )

# ----------------------------------------------------------------------- #
# 03 - MD Simulations - 5ns Production with 4fs timestep
# ----------------------------------------------------------------------- #
# Simulation Object
integrator = openmm.LangevinMiddleIntegrator(
    temperature, 1.0 / unit.picoseconds, 2.0 * unit.femtoseconds
)
simulation = app.Simulation(
    prmtop.topology,
    system,
    integrator,
    openmm.Platform.getPlatformByName("OpenCL"),
)
simulation.loadState(f"{work_dir}/checkpoint_equilibration.xml")

# Reporters
simulation.reporters.append(
    app.CheckpointReporter(f"{work_dir}/checkpoint_production.xml", 25000, writeState=True)
)
simulation.reporters.append(
    app.DCDReporter(f"{work_dir}/trajectory_production.dcd", 2500)
)
simulation.reporters.append(
    app.StateDataReporter(
        f"{work_dir}/reporter_production.csv",
        2500,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
    )
)

# MD Simulation
simulation.step(2500000)

# Write final snapshot
simulation.saveState(f"{work_dir}/checkpoint_production.xml")
with open(f"{work_dir}/complex_production.pdb", "w") as f:
    app.PDBFile.writeFile(
        simulation.topology,
        simulation.context.getState(getPositions=True).getPositions(),
        f,
        keepIds=True,
    )

