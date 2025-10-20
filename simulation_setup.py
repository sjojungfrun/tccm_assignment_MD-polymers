"""
OpenMM Simulation Setup Module

This module contains functions for setting up molecular dynamics simulations,
including force field creation, system building, and initial configuration.
"""

from openmm.app import *
from openmm import *
from openmm.unit import *
from openmm import unit
import numpy as np


def create_forcefield(ff_files=['amber14-all.xml', 'amber14/tip3p.xml']):
    """
    Create a force field object with specified force field files.
    
    Parameters:
    -----------
    ff_files : list of str
        List of force field XML files to load
        
    Returns:
    --------
    ForceField
        OpenMM ForceField object
    """
    return ForceField(*ff_files)


def estimate_box_size(n_waters, volume_per_water=30):
    """
    Estimate cubic box size for a given number of water molecules.
    
    Parameters:
    -----------
    n_waters : int
        Number of water molecules
    volume_per_water : float
        Volume per water molecule in Ų (default: 30 Ų for liquid water)
        
    Returns:
    --------
    Quantity
        Box length with units
    """
    total_volume = n_waters * volume_per_water
    box_length_angstroms = total_volume ** (1/3)
    return box_length_angstroms * unit.angstroms


def create_water_box(forcefield, n_waters, box_length=None, water_model='tip3p'):
    """
    Create a water box with specified number of molecules or box size.
    
    Parameters:
    -----------
    forcefield : ForceField
        OpenMM ForceField object
    n_waters : int
        Target number of water molecules (used for box size estimation if box_length is None)
    box_length : Quantity, optional
        Box length with units. If None, estimated from n_waters
    water_model : str
        Water model to use (default: 'tip3p')
        
    Returns:
    --------
    tuple
        (modeller, actual_box_length) where modeller is the Modeller object
        and actual_box_length is the box length with units
    """
    # Create empty modeller
    modeller = Modeller(Topology(), [])
    
    # Estimate box size if not provided
    if box_length is None:
        box_length = estimate_box_size(n_waters)
    
    # Add solvent in a cubic box
    box_size = Vec3(
        box_length.value_in_unit(unit.angstroms),
        box_length.value_in_unit(unit.angstroms),
        box_length.value_in_unit(unit.angstroms)
    ) * unit.angstroms
    
    modeller.addSolvent(forcefield, model=water_model, boxSize=box_size)
    
    return modeller, box_length


def create_system(forcefield, topology, nonbonded_method=PME, 
                 cutoff=0.5*unit.nanometers, constraints=HBonds,
                 temperature=300*unit.kelvin, pressure=1*unit.bar):
    """
    Create an OpenMM system with specified parameters.
    
    Parameters:
    -----------
    forcefield : ForceField
        OpenMM ForceField object
    topology : Topology
        System topology
    nonbonded_method : method
        Nonbonded interaction method (default: PME)
    cutoff : Quantity
        Nonbonded cutoff distance
    constraints : constraint type
        Bond constraints (default: HBonds)
    temperature : Quantity
        Temperature for barostat
    pressure : Quantity
        Pressure for barostat
        
    Returns:
    --------
    System
        OpenMM System object with barostat added
    """
    # Create system
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=cutoff,
        constraints=constraints
    )
    
    return system


def create_integrator(temperature=300*unit.kelvin, friction=1/unit.picoseconds, 
                     timestep=2*unit.femtoseconds):
    """
    Create a Langevin integrator for molecular dynamics.
    
    Parameters:
    -----------
    temperature : Quantity
        Target temperature
    friction : Quantity
        Friction coefficient
    timestep : Quantity
        Integration timestep
        
    Returns:
    --------
    LangevinMiddleIntegrator
        OpenMM integrator object
    """
    return LangevinMiddleIntegrator(temperature, friction, timestep)


def setup_simulation(n_waters=50, ensemble='nvt', temperature=300*unit.kelvin, 
                    pressure=1*unit.bar, timestep=2*unit.femtoseconds, 
                    cutoff=0.5*unit.nanometers,
                    ff_files=['amber14-all.xml', 'amber14/tip3p.xml']):
    """
    Complete simulation setup function that creates all necessary components.
    
    Parameters:
    -----------
    n_waters : int
        Target number of water molecules
    ensemble : string
        Ensemble of the simulation, to choose between 'nvt' or 'npt'
    temperature : Quantity
        Simulation temperature
    pressure : Quantity
        Simulation pressure
    timestep : Quantity
        Integration timestep
    cutoff : Quantity
        Nonbonded cutoff distance
    ff_files : list
        Force field XML files
        
    Returns:
    --------
    dict
        Dictionary containing 'simulation', 'modeller', 'box_length', and 'system'
    """
    print(f"Setting up water box with {n_waters} molecules...")
    
    # Create force field
    forcefield = create_forcefield(ff_files)
    
    # Create water box
    modeller, box_length = create_water_box(forcefield, n_waters)
    
    print(f"Box dimensions: {box_length.value_in_unit(unit.nanometers):.2f} nm x "
          f"{box_length.value_in_unit(unit.nanometers):.2f} nm x "
          f"{box_length.value_in_unit(unit.nanometers):.2f} nm")
    print(f"Actual number of water molecules: {modeller.topology.getNumResidues()}")
    
    # Create system
    system = create_system(forcefield, modeller.topology, cutoff=cutoff,
                          temperature=temperature, pressure=pressure)
    
    # Create integrator
    integrator = create_integrator(temperature=temperature, timestep=timestep)
    
    # Add barostat if NPT ensemble
    ensemble = ensemble.lower()
    if ensemble == 'npt':
        system.addForce(MonteCarloBarostat(pressure, temperature))
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    return {
        'simulation': simulation,
        'modeller': modeller,
        'box_length': box_length,
        'system': system
    }