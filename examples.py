"""
Example usage of the modular OpenMM simulation library

This script demonstrates how to use the modular functions to create
custom simulation protocols with different parameters.
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit

from simulation_setup import setup_simulation, create_forcefield, create_water_box, create_system
from simulation_runner import run_complete_simulation, run_custom_protocol
from analysis import complete_analysis

def example_basic_simulation():
    """Example of a basic water density simulation."""
    print("Running basic water density simulation...")
    
    # Setup with different parameters
    simulation_components = setup_simulation(
        n_waters=30,  # Smaller system
        temperature=298*kelvin,  # Room temperature
        pressure=1*bar,
        timestep=1*unit.femtoseconds,  # Smaller timestep
        cutoff=0.4*unit.nanometers
    )
    
    # Run simulation
    results = run_complete_simulation(
        simulation=simulation_components['simulation'],
        nvt_steps=2000,
        npt_steps=3000,
        production_steps=10000,
        temperature=298*kelvin,
        output_log='example_basic.log',
        output_pdb='example_basic.pdb'
    )
    
    # Analyze
    n_molecules = simulation_components['modeller'].topology.getNumResidues()
    analysis = complete_analysis(
        volumes=results['volumes'],
        n_molecules=n_molecules,
        output_file='example_basic_analysis.txt'
    )
    
    return analysis

def example_custom_protocol():
    """Example of a custom simulation protocol."""
    print("\nRunning custom protocol simulation...")
    
    # Setup simulation
    simulation_components = setup_simulation(n_waters=25)
    simulation = simulation_components['simulation']
    
    # Define custom protocol
    protocol = [
        {'type': 'minimize', 'max_iterations': 200},
        {'type': 'equilibrate', 'steps': 1000},  # Short NVT
        {'type': 'equilibrate', 'steps': 2000},  # Short NPT  
        {'type': 'production', 'steps': 5000, 'collect_data': True, 'data_interval': 50}
    ]
    
    # Run custom protocol
    results = run_custom_protocol(simulation, protocol)
    
    # Analyze results
    if results['volumes']:
        n_molecules = simulation_components['modeller'].topology.getNumResidues()
        analysis = complete_analysis(
            volumes=results['volumes'],
            n_molecules=n_molecules,
            output_file='example_custom_analysis.txt'
        )
        return analysis
    
    return None

def example_step_by_step():
    """Example of building simulation step by step."""
    print("\nRunning step-by-step simulation setup...")
    
    # Step 1: Create force field
    ff = create_forcefield(['amber14-all.xml', 'amber14/tip3p.xml'])
    
    # Step 2: Create water box
    modeller, box_length = create_water_box(ff, n_waters=20)
    print(f"Created box with {modeller.topology.getNumResidues()} molecules")
    
    # Step 3: Create system
    system = create_system(ff, modeller.topology, cutoff=0.4*unit.nanometers)
    
    # Step 4: Create integrator and simulation
    from simulation_setup import create_integrator
    integrator = create_integrator(temperature=295*kelvin)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Step 5: Run with custom parameters
    results = run_complete_simulation(
        simulation=simulation,
        nvt_steps=1000,
        npt_steps=1500,
        production_steps=5000,
        temperature=295*kelvin,
        output_log='example_stepby_step.log',
        output_pdb='example_step_by_step.pdb'
    )
    
    return results

if __name__ == "__main__":
    print("OpenMM Modular Simulation Examples")
    print("="*40)
    
    # Run examples
    try:
        basic_results = example_basic_simulation()
        print(f"Basic simulation completed. Density: {basic_results['density_results']['density']:.4f} g/cm³")
    except Exception as e:
        print(f"Basic simulation failed: {e}")
    
    try:
        custom_results = example_custom_protocol()
        if custom_results:
            print(f"Custom protocol completed. Density: {custom_results['density_results']['density']:.4f} g/cm³")
    except Exception as e:
        print(f"Custom protocol failed: {e}")
    
    try:
        step_results = example_step_by_step()
        print(f"Step-by-step simulation completed. {len(step_results['volumes'])} volume samples collected.")
    except Exception as e:
        print(f"Step-by-step simulation failed: {e}")
    
    print("\nAll examples completed!")