import openmm
import openmm.app as app
import openmm.unit as unit

from hard_sphere import make_guest_wca_repulsive

# User Option
# ----------------------------------------------------------------------- #
guest_resname = "DM1"
host_resname = "CB8"
solvent_resname = "HOH"

hard_sphere_sigma = 3.0 * unit.angstrom
hard_sphere_epsilon = 0.1 * unit.kilocalorie_per_mole
hard_sphere_radius = 5.0 * unit.angstrom

wca_repulsive = 50
wca_attractive = 49

# Load Amber Files
# ----------------------------------------------------------------------- #
prmtop = app.AmberPrmtopFile("structures/CB8-Guest-sol.prmtop")
inpcrd = app.AmberInpcrdFile("structures/CB8-Guest-sol.rst7")

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

# Set Host molecule either "hard-sphere" or weakly dispersive
# ----------------------------------------------------------------------- #
# make_host_wca_dispersive(
#    system,
#    host,
#    solvent,
#    lambda_value=0.0,
#    force_group=9,
# )

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

# Set Host charges and LJ parameters to zero for weakly dispersive host
# ----------------------------------------------------------------------- #
# for atom in host:
#    charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
#    nonbonded.setParticleParameters(atom, 0.0, sigma, 0.0)
#
# for i_exception in range(nonbonded.getNumExceptions()):
#    atom_i, atom_j, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
#        i_exception
#    )
#    if atom_i in host and atom_j in host:
#        nonbonded.setExceptionParameters(i_exception, atom_i, atom_j, 0.0, sigma, 0.0)

with open("system_hard-sphere.xml", "w") as file:
    file.write(openmm.XmlSerializer.serialize(system))
