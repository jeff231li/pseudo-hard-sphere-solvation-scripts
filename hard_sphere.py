import openmm
import openmm.app as app
import openmm.unit as unit


def make_guest_wca_repulsive(
    system,
    guest,
    host,
    solvents,
    particle_radius=5.0 * unit.angstrom,
    wca_sigma=3.0 * unit.angstrom,
    wca_epsilon=1.0 * unit.kilocalorie_per_mole,
    coeff_r=50,
    coeff_a=49,
    force_group=10,
):
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    ][0]

    wca_repulsive = openmm.CustomNonbondedForce(
        "U_repulsive;"
        "U_repulsive = step(R_particle - r) * (U_Mie + epsilon_wall);"
        "U_Mie = prefactor * epsilon_wall * ((sigma_wall/r_prime)^coeff_r - (sigma_wall/r_prime)^coeff_a);"
        "prefactor = coeff_r/(coeff_r-coeff_a) * (coeff_r/coeff_a)^(coeff_r/(coeff_r-coeff_a));"
        "r_prime = r - (R_particle - R_min);"
        "R_min = sigma_wall * (coeff_r/coeff_a)^(1/(coeff_r-coeff_a));"
    )
    wca_repulsive.addGlobalParameter("R_particle", particle_radius)
    wca_repulsive.addGlobalParameter("coeff_r", coeff_r)
    wca_repulsive.addGlobalParameter("coeff_a", coeff_a)
    wca_repulsive.addGlobalParameter("sigma_wall", wca_sigma)
    wca_repulsive.addGlobalParameter("epsilon_wall", wca_epsilon)
    wca_repulsive.addPerParticleParameter("sigma")
    wca_repulsive.addPerParticleParameter("epsilon")
    wca_repulsive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    wca_repulsive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_repulsive.setUseLongRangeCorrection(False)
    wca_repulsive.setForceGroup(force_group)

    wca_repulsive.addInteractionGroup(guest, (host + solvents))

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        if atom in guest:
            wca_repulsive.addParticle([wca_sigma, wca_epsilon])
        else:
            wca_repulsive.addParticle([sigma, epsilon])

    system.addForce(wca_repulsive)

    # Set LJ parameters to zero
    for atom in guest:
        [charge, sigma, epsilon] = nonbonded.getParticleParameters(atom)
        nonbonded.setParticleParameters(atom, 0.0, wca_sigma, 0.0)

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
            exception_index
        )
        wca_repulsive.addExclusion(iatom, jatom)


def make_host_wca_dispersive(
    system,
    host_atoms,
    solvent_atoms,
    lambda_scale=1,
    force_group=10,
):
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    ][0]

    energy_expression = (
        "U_dispersive;"
        "U_dispersive = step(R_min - r)*(U_LJ + epsilon*(1 - lambda_dispersive)) + step(r - R_min)*(lambda_dispersive * U_LJ);"
        "U_LJ = 4 * epsilon * x * (x - 1.0);"
        "x = (sigma / r)^6;"
        "R_min = sigma * 2^(1/6);"
        "sigma = 0.5*(sigma1 + sigma2);"
        "epsilon = sqrt(epsilon1 * epsilon2);"
    )

    wca_dispersive = openmm.CustomNonbondedForce(energy_expression)
    wca_dispersive.addGlobalParameter("lambda_dispersive", lambda_scale)
    wca_dispersive.addPerParticleParameter("sigma")
    wca_dispersive.addPerParticleParameter("epsilon")
    wca_dispersive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    wca_dispersive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_dispersive.setUseLongRangeCorrection(nonbonded.getUseDispersionCorrection())
    wca_dispersive.addInteractionGroup(host_atoms, solvent_atoms)
    wca_dispersive.setForceGroup(force_group)

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        wca_dispersive.addParticle([sigma, epsilon])

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
            exception_index
        )
        wca_dispersive.addExclusion(iatom, jatom)

    # Add force to System
    system.addForce(wca_dispersive)


def alchemicalize_guest_wca_repulsive(
    system,
    guest,
    solvents,
    particle_radius=5.0 * unit.angstrom,
    wca_sigma=3.0 * unit.angstrom,
    wca_epsilon=1.0 * unit.kilocalorie_per_mole,
    coeff_r=50,
    coeff_a=49,
    softcore_alpha=0.5,
    lambda_value=1.0,
    force_group=10,
):
    nonbonded = [
        force
        for force in system.getForces()
        if isinstance(force, openmm.NonbondedForce)
    ][0]

    wca_repulsive = openmm.CustomNonbondedForce(
        "lambda_sterics * U_repulsive;"
        "U_repulsive = step(R_particle - r) * (U_Mie + epsilon_wall);"
        "U_Mie = prefactor * epsilon_wall * (1/repulsive - 1/dispersive);"
        "repulsive = (dispersive)^(coeff_r/coeff_a);"
        "dispersive = softcore_alpha*(1.0-lambda_sterics) + (r_prime/sigma_wall)^coeff_a;"
        "prefactor = (coeff_r/(coeff_r - coeff_a)) * (coeff_r/coeff_a)^(coeff_a/(coeff_r-coeff_a));"
        "r_prime = r - (R_particle - R_min);"
        "R_min = sigma_wall * ((coeff_r/coeff_a)^(coeff_a/(coeff_r-coeff_a)) - softcore_alpha*(1.0-lambda_sterics))^(1/coeff_a);"
    )
    wca_repulsive.addGlobalParameter("lambda_sterics", lambda_value)
    wca_repulsive.addGlobalParameter("softcore_alpha", softcore_alpha)
    wca_repulsive.addGlobalParameter("R_particle", particle_radius)
    wca_repulsive.addGlobalParameter("coeff_r", coeff_r)
    wca_repulsive.addGlobalParameter("coeff_a", coeff_a)
    wca_repulsive.addGlobalParameter("sigma_wall", wca_sigma)
    wca_repulsive.addGlobalParameter("epsilon_wall", wca_epsilon)
    wca_repulsive.addPerParticleParameter("sigma")
    wca_repulsive.addPerParticleParameter("epsilon")
    wca_repulsive.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
    wca_repulsive.setCutoffDistance(nonbonded.getCutoffDistance())
    wca_repulsive.setUseLongRangeCorrection(False)
    wca_repulsive.setForceGroup(force_group)

    wca_repulsive.addInteractionGroup(guest, solvents)

    # Set LJ parameters
    for atom in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(atom)
        if atom in guest:
            wca_repulsive.addParticle([wca_sigma, wca_epsilon])
        else:
            wca_repulsive.addParticle([sigma, epsilon])

    # Set LJ parameters to zero
    for atom in guest:
        [charge, sigma, epsilon] = nonbonded.getParticleParameters(atom)
        nonbonded.setParticleParameters(atom, 0.0, wca_sigma, 0.0)

    # Transfer Exclusion
    for exception_index in range(nonbonded.getNumExceptions()):
        iatom, jatom, chargeprod, sigma, epsilon = nonbonded.getExceptionParameters(
            exception_index
        )
        wca_repulsive.addExclusion(iatom, jatom)

    # Add Force to System
    system.addForce(wca_repulsive)
