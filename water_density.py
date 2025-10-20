"""
OpenMM script to create a water box, equilibrate it, and calculate density
Using AMBER force field (TIP3P water model) - Modular Version
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit
import numpy as np

# Import our custom modules
from simulation_setup import setup_simulation
from simulation_runner import run_complete_simulation
from analysis import complete_analysis

# Parameters
n_waters = 50
temperature = 300 * kelvin
pressure = 1 * bar
nvt_equilibration_steps = 5000
npt_equilibration_steps = 10000
production_steps = 50000
timestep = 2 * unit.femtoseconds
cutoff = 0.5 * unit.nanometers

def main():
    """Main function to run the water density simulation."""
    
    # Setup simulation
    simulation_components = setup_simulation(
        n_waters=n_waters,
        temperature=temperature,
        pressure=pressure,
        timestep=timestep,
        cutoff=cutoff
    )
    
    simulation = simulation_components['simulation']
    modeller = simulation_components['modeller']
    
    # Run complete simulation
    simulation_results = run_complete_simulation(
        simulation=simulation,
        nvt_steps=nvt_equilibration_steps,
        npt_steps=npt_equilibration_steps,
        production_steps=production_steps,
        temperature=temperature,
        data_collection_interval=100,
        output_log='water_output.log',
        output_pdb='final_water_box.pdb'
    )
    
    # Analyze results
    n_molecules = modeller.topology.getNumResidues()
    volumes = simulation_results['volumes']
    
    analysis_results = complete_analysis(
        volumes=volumes,
        n_molecules=n_molecules,
        reference_density=1.0,
        save_to_file=True,
        output_file='analysis_results.txt'
    )
    
    return analysis_results

if __name__ == "__main__":
    results = main()