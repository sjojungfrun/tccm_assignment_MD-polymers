import os

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit

# Import our custom modules
from simulation_setup import setup_simulation
from simulation_runner import minimize_energy, equilibrate_system, run_production_vol, save_final_structure
from analysis import complete_analysis


def _temperature_suffix(temperature_quantity):
    """Return suffix like '300K' for a given temperature."""
    if hasattr(temperature_quantity, 'value_in_unit'):
        kelvin_value = temperature_quantity.value_in_unit(unit.kelvin)
    else:
        kelvin_value = temperature_quantity
    return f"{int(kelvin_value):d}K"


def _with_temperature_suffix(filename, temperature_quantity):
    """Append the temperature suffix before a filename extension."""
    name, ext = os.path.splitext(filename)
    return f"{name}_{_temperature_suffix(temperature_quantity)}{ext}"


# Parameters

def density_at_temperature(
    n_waters = 50,
    temperature = 300 * kelvin,
    pressure = 1 * bar,
    npt_equilibration_steps = 1000,
    production_steps = 50000,
    timestep = 2 * unit.femtoseconds,
    cutoff = 0.5 * unit.nanometers,
    reference_density = 1.0  # g/cm³ for water at room temperature
):
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
    
    # Perform first an energy minimization
    minimize_energy(simulation, max_iterations=500)

    # Perform first system equiliration on NVT ensemble
    print("Starting NPT equilibration...")
    equilibrate_system(
        simulation,
        total_steps=npt_equilibration_steps,
        ensemble='npt',
        temperature=temperature,
        pressure=pressure
    )

    output_log = _with_temperature_suffix('water_output.log', temperature)
    output_pdb = _with_temperature_suffix('final_water_box.pdb', temperature)
    analysis_output_file = _with_temperature_suffix('analysis_results.txt', temperature)
    analysis_json_file = _with_temperature_suffix('analysis_results.json', temperature)
    summary_filename = _with_temperature_suffix('density_summary.txt', temperature)

    # Perform production run in NPT ensemble and collect volume data
    sim_results = run_production_vol(
        simulation,
        ensemble='npt',
        temperature=temperature,
        pressure=pressure,
        production_steps=production_steps,
        data_collection_interval=100,
        output_file=output_log
    )
    # Save final structure
    save_final_structure(
        simulation,
        output_file=output_pdb
    )
    
    # Analyze results
    n_molecules = modeller.topology.getNumResidues()
    volumes = sim_results['volumes']
    
    analysis_results = complete_analysis(
        volumes=volumes,
        n_molecules=n_molecules,
        reference_density=reference_density,
        save_to_file=True,
        output_file=analysis_output_file,
        json_file=analysis_json_file
    )

    density_value = analysis_results['density_results']['density']
    rel_error = analysis_results['relative_error']
    compressibility = analysis_results['compressibility']

    print(f"Calculated density: {density_value:.4f} g/cm³")
    print(f"Relative error vs reference: {rel_error:.2f}%")
    print(f"Estimated compressibility: {compressibility:.3e} 1/Pa")

    summary_lines = [
        f"Average volume: {analysis_results['density_results']['avg_volume_nm3']:.4f} ± "
        f"{analysis_results['density_results']['std_volume_nm3']:.4f} nm³",
        f"Average density: {density_value:.4f} g/cm³",
        f"Expected density (experimental): {reference_density:.1f} g/cm³",
        f"Relative error: {rel_error:.2f}%",
        "=" * 50,
    ]

    with open(summary_filename, 'w', encoding='utf-8') as summary_file:
        summary_file.write("\n".join(summary_lines))
    print(f"Summary written to '{summary_filename}'")
    print(f"Simulation log: {output_log}")
    print(f"Final structure: {output_pdb}")
    print(f"Analysis files: {analysis_output_file}, {analysis_json_file}")
    
    return analysis_results

if __name__ == "__main__":
    density_at_temperature(
        n_waters = 50,
        temperature = 300 * kelvin,
        pressure = 1 * bar,
        npt_equilibration_steps = 1000,
        production_steps = 50000,
        timestep = 2 * unit.femtoseconds,
        cutoff = 0.5 * unit.nanometers,
        reference_density = 1.0  # g/cm³ for water at room temperature
    )
