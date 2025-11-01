from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit
import numpy as np


def minimize_energy(simulation, max_iterations=100, tolerance=10*unit.kilojoules_per_mole/unit.nanometers):
    """
    Minimize the energy of the system.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    max_iterations : int
        Maximum number of minimization iterations
    tolerance : Quantity
        Energy tolerance for convergence
    """
    print("\nMinimizing energy...")
    simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)


def equilibrate_system(simulation, total_steps, ensemble='npt', 
                      temperature=300*unit.kelvin, pressure=1*unit.bar):
    """
    Equilibrate the system in the specified ensemble.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    total_steps : int
        Number of equilibration steps
    ensemble : str
        Either 'nvt' or 'npt'
    temperature : Quantity
        Target temperature
    pressure : Quantity
        Target pressure (only used for NPT)
    """
    ensemble = ensemble.lower()
    print(f"Equilibrating in {ensemble.upper()} ensemble for {total_steps} steps...")
    
    # Set velocities from temperature
    simulation.context.setVelocitiesToTemperature(temperature)

    # Set simulation to proper ensemble
    switch_ensemble(simulation, ensemble, temperature, pressure)

    # Run equilibration
    simulation.step(total_steps)
    print(f"Equilibration complete: {total_steps} steps")
    
def switch_ensemble(simulation, ensemble='npt', temperature=300*unit.kelvin, pressure=1*unit.bar):
    # Handle barostat for NPT
    system = simulation.system
    barostat_present = False
    barostat_index = None
    
    # Check if barostat exists
    for index in range(system.getNumForces()):
        force = system.getForce(index)
        if isinstance(force, MonteCarloBarostat):
            barostat_present = True
            barostat_index = index
            break
    
    if ensemble == 'npt':
        if not barostat_present:
            # Add barostat for NPT
            system.addForce(MonteCarloBarostat(pressure, temperature))
            simulation.context.reinitialize(preserveState=True)
            print(f"Added MonteCarloBarostat at {pressure} and {temperature}")
    elif ensemble == 'nvt':
        if barostat_present:
            # Remove barostat for NVT
            system.removeForce(barostat_index)
            simulation.context.reinitialize(preserveState=True)
            print("Removed MonteCarloBarostat for NVT ensemble")
    else:
        raise ValueError(f"Invalid ensemble '{ensemble}'. Must be 'nvt' or 'npt'.")
    


def run_production_vol(simulation, ensemble='npt',  
                      temperature=300*unit.kelvin, pressure=1*unit.bar,
                      production_steps=50000, data_collection_interval=100,
                      output_file='simulation_output.log', log_interval=1000):
    """
    Run production simulation and collect volume data.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    production_steps : int
        Number of production steps
    data_collection_interval : int
        Interval for collecting volume data
    output_file : str
        Name of the output log file
    log_interval : int
        Interval for writing log data
        
    Returns:
    --------
    list
        List of volumes in nm³
    """
    print(f"\nRunning production for {production_steps} steps...")
    
    # Add reporter for logging
    simulation.reporters.append(StateDataReporter(output_file, log_interval,
        step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
    
    volumes = []
    
    # Set simulation to proper ensemble
    switch_ensemble(simulation, ensemble, temperature, pressure)

    # Collect volume data during production
    for _ in range(production_steps // data_collection_interval):
        simulation.step(data_collection_interval)
        state = simulation.context.getState(enforcePeriodicBox=True)
        box_vectors = state.getPeriodicBoxVectors()
        
        # Calculate volume from box vectors
        a, b, c = box_vectors
        volume = a[0] * b[1] * c[2]  # For cubic/orthorhombic box
        volumes.append(volume.value_in_unit(unit.nanometers**3))
    
    print(f"Production complete: {production_steps} steps")
    print(f"Data saved to '{output_file}'")

    return {
        'volumes': volumes,
        'production_steps': production_steps,
        'ensemble': ensemble
    }


def save_final_structure(simulation, output_file='final_structure.pdb'):
    """
    Save the final configuration to a PDB file.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    output_file : str
        Name of the output PDB file
        
    Returns:
    --------
    str
        Path to the saved file
    """
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(output_file, 'w'))
    print(f"\nFinal configuration saved to '{output_file}'")
    return output_file


def run_complete_simulation(simulation, equilibration_steps=10000, 
                           production_steps=50000, ensemble='npt',
                           temperature=300*unit.kelvin, pressure=1*unit.bar,
                           data_collection_interval=100,
                           output_log='simulation_output.log', 
                           output_pdb='final_structure.pdb'):
    """
    Run a complete simulation including minimization, equilibration, and production.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    equilibration_steps : int
        Number of equilibration steps
    production_steps : int
        Number of production steps
    ensemble : str
        Either 'nvt' or 'npt'
    temperature : Quantity
        Target temperature
    pressure : Quantity
        Target pressure (only used for NPT)
    data_collection_interval : int
        Interval for collecting volume data
    output_log : str
        Name of the output log file
    output_pdb : str
        Name of the output PDB file
        
    Returns:
    --------
    dict
        Dictionary containing 'volumes', 'equilibration_steps', 'production_steps'
    """
    # Energy minimization
    minimize_energy(simulation)
    
    # Equilibration in specified ensemble
    equilibrate_system(simulation, equilibration_steps, ensemble, temperature, pressure)
    
    # Production run
    volumes = run_production(simulation, production_steps, data_collection_interval, output_log)
    
    # Save final structure
    save_final_structure(simulation, output_pdb)
    
    print(f"\nSimulation complete!")
    print(f"Data saved to '{output_log}'")
    print(f"Structure saved to '{output_pdb}'")
    
    return {
        'volumes': volumes,
        'equilibration_steps': equilibration_steps,
        'production_steps': production_steps,
        'ensemble': ensemble
    }