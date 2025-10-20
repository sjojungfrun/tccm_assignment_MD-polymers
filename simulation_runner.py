"""
OpenMM Simulation Runner Module

This module contains functions for running molecular dynamics simulations,
including energy minimization, equilibration, production runs, and data collection.
"""

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
        
    Returns:
    --------
    None
    """
    print("\nMinimizing energy...")
    simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)


def equilibrate_system(simulation, total_steps=5000, temperature=300*unit.kelvin):
    """
    Equilibrate the system with NVT followed by NPT.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    total_steps : int
        Number of steps of simulation
    temperature : Quantity
        Target temperature
        
    Returns:
    --------
    int
        Total number of equilibration steps performed
    """
    print(f"Equilibrating for {total_steps} steps...")
    
    # Set initial velocities
    simulation.context.setVelocitiesToTemperature(temperature)
    
    # Run equilibration
    simulation.step(total_steps)
    
    return total_steps


def run_production(simulation, production_steps=50000, data_collection_interval=100,
                  output_file='water_output.log', log_interval=1000):
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
    print(f"Running production for {production_steps} steps...")
    
    # Add reporter for logging
    simulation.reporters.append(StateDataReporter(output_file, log_interval,
        step=True, potentialEnergy=True, temperature=True, volume=True, density=True))
    
    volumes = []
    
    # Collect volume data during production
    for _ in range(production_steps // data_collection_interval):
        simulation.step(data_collection_interval)
        state = simulation.context.getState()
        box_vectors = state.getPeriodicBoxVectors()
        
        # Calculate volume from box vectors
        a, b, c = box_vectors
        volume = a[0] * b[1] * c[2]  # For cubic/orthorhombic box
        volumes.append(volume.value_in_unit(unit.nanometers**3))
    
    return volumes


def save_final_structure(simulation, output_file='final_water_box.pdb'):
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


def run_complete_simulation(simulation, total_steps=5000, prod_steps=50000, 
                           temperature=300*unit.kelvin, data_collection_interval=100, 
                           output_log='water_output.log', output_pdb='final_water_box.pdb'):
    """
    Run a complete simulation including minimization, equilibration, and production.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    eq_steps : int
        Number of total equilibration steps
    prod_steps : int
        Number of production steps
    temperature : Quantity
        Target temperature
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
    
    # Equilibration
    total_equilibration_steps = equilibrate_system(simulation, total_steps, temperature)
    
    # Production run
    volumes = run_production(simulation, prod_steps, data_collection_interval, output_log)
    
    # Save final structure
    save_final_structure(simulation, output_pdb)
    
    print(f"Simulation data saved to '{output_log}'")
    
    return {
        'volumes': volumes,
        'equilibration_steps': total_equilibration_steps,
        'production_steps': prod_steps
    }


def run_custom_protocol(simulation, protocol_steps, temperature=300*unit.kelvin):
    """
    Run a custom simulation protocol with specified steps.
    
    Parameters:
    -----------
    simulation : Simulation
        OpenMM Simulation object
    protocol_steps : list of dict
        List of protocol steps, each containing 'type', 'steps', and optional parameters
        Example: [{'type': 'minimize', 'max_iterations': 100},
                 {'type': 'equilibrate', 'steps': 1000},
                 {'type': 'production', 'steps': 5000, 'collect_data': True}]
    temperature : Quantity
        Target temperature
        
    Returns:
    --------
    dict
        Results from the protocol execution
    """
    results = {'volumes': [], 'total_steps': 0}
    
    for i, step in enumerate(protocol_steps):
        step_type = step['type'].lower()
        
        if step_type == 'minimize':
            max_iter = step.get('max_iterations', 100)
            minimize_energy(simulation, max_iter)
            
        elif step_type == 'equilibrate':
            steps = step['steps']
            if i == 0:  # First equilibration step
                simulation.context.setVelocitiesToTemperature(temperature)
            simulation.step(steps)
            results['total_steps'] += steps
            print(f"Equilibration step {i+1}: {steps} steps completed")
            
        elif step_type == 'production':
            steps = step['steps']
            collect_data = step.get('collect_data', False)
            interval = step.get('data_interval', 100)
            
            if collect_data:
                volumes = run_production(simulation, steps, interval)
                results['volumes'].extend(volumes)
            else:
                simulation.step(steps)
            
            results['total_steps'] += steps
            print(f"Production step {i+1}: {steps} steps completed")
    
    return results