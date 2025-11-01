"""
OpenMM Analysis Module

This module contains functions for analyzing simulation results,
including density calculations, statistical analysis, and result reporting.
"""

import numpy as np
from openmm.unit import *
from openmm import unit


def calculate_water_density(volumes, n_molecules, mass_per_molecule=18.015):
    """
    Calculate water density from volume data.
    
    Parameters:
    -----------
    volumes : list
        List of volumes in nm³
    n_molecules : int
        Number of water molecules
    mass_per_molecule : float
        Mass per water molecule in g/mol (default: 18.015)
        
    Returns:
    --------
    dict
        Dictionary containing density statistics and volume statistics
    """
    # Calculate volume statistics
    avg_volume_nm3 = np.mean(volumes)
    std_volume_nm3 = np.std(volumes)
    
    # Calculate mass
    avogadro = 6.022e23  # molecules/mol
    total_mass = n_molecules * mass_per_molecule  # g/mol
    total_mass_g = total_mass / avogadro  # g
    
    # Convert volume to cm³ (1 nm³ = 1e-21 cm³)
    avg_volume_cm3 = avg_volume_nm3 * 1e-21
    
    # Calculate density in g/cm³
    density = total_mass_g / avg_volume_cm3
    
    return {
        'density': density,
        'avg_volume_nm3': avg_volume_nm3,
        'std_volume_nm3': std_volume_nm3,
        'avg_volume_cm3': avg_volume_cm3,
        'total_mass_g': total_mass_g,
        'n_molecules': n_molecules
    }


def calculate_density_error(calculated_density, reference_density=1.0):
    """
    Calculate relative error in density.
    
    Parameters:
    -----------
    calculated_density : float
        Calculated density in g/cm³
    reference_density : float
        Reference density in g/cm³ (default: 1.0 for water)
        
    Returns:
    --------
    float
        Relative error as percentage
    """
    return abs(calculated_density - reference_density) / reference_density * 100


def analyze_volume_fluctuations(volumes):
    """
    Analyze volume fluctuations during simulation.
    
    Parameters:
    -----------
    volumes : list
        List of volumes in nm³
        
    Returns:
    --------
    dict
        Dictionary containing fluctuation statistics
    """
    volumes_array = np.array(volumes)
    
    return {
        'mean': np.mean(volumes_array),
        'std': np.std(volumes_array),
        'min': np.min(volumes_array),
        'max': np.max(volumes_array),
        'range': np.max(volumes_array) - np.min(volumes_array),
        'cv': np.std(volumes_array) / np.mean(volumes_array) * 100,  # Coefficient of variation
        'n_samples': len(volumes_array)
    }


def calculate_compressibility(volumes, temperature=300, kb=1.381e-23):
    """
    Estimate isothermal compressibility from volume fluctuations.
    
    Parameters:
    -----------
    volumes : list
        List of volumes in nm³
    temperature : float
        Temperature in Kelvin
    kb : float
        Boltzmann constant in J/K
        
    Returns:
    --------
    float
        Isothermal compressibility in 1/Pa
    """
    volumes_array = np.array(volumes) * 1e-27  # Convert nm³ to m³
    avg_volume = np.mean(volumes_array)
    volume_variance = np.var(volumes_array)
    
    # κT = <ΔV²> / (kB * T * <V>)
    compressibility = volume_variance / (kb * temperature * avg_volume)
    
    return compressibility


def print_analysis_results(analysis_results, reference_density=1.0):
    """
    Print formatted analysis results.
    
    Parameters:
    -----------
    analysis_results : dict
        Results from calculate_water_density function
    reference_density : float
        Reference density for comparison
        
    Returns:
    --------
    None
    """
    density = analysis_results['density']
    error = calculate_density_error(density, reference_density)
    
    print("\n" + "="*50)
    print("ANALYSIS RESULTS")
    print("="*50)
    print(f"Number of water molecules: {analysis_results['n_molecules']}")
    print(f"Average volume: {analysis_results['avg_volume_nm3']:.4f} ± {analysis_results['std_volume_nm3']:.4f} nm³")
    print(f"Average density: {density:.4f} g/cm³")
    print(f"Expected density (experimental): {reference_density:.1f} g/cm³")
    print(f"Relative error: {error:.2f}%")
    print("="*50)


def print_volume_statistics(volume_stats):
    """
    Print detailed volume fluctuation statistics.
    
    Parameters:
    -----------
    volume_stats : dict
        Results from analyze_volume_fluctuations function
        
    Returns:
    --------
    None
    """
    print("\nVOLUME FLUCTUATION ANALYSIS")
    print("-" * 30)
    print(f"Mean volume: {volume_stats['mean']:.4f} nm³")
    print(f"Standard deviation: {volume_stats['std']:.4f} nm³")
    print(f"Min volume: {volume_stats['min']:.4f} nm³")
    print(f"Max volume: {volume_stats['max']:.4f} nm³")
    print(f"Volume range: {volume_stats['range']:.4f} nm³")
    print(f"Coefficient of variation: {volume_stats['cv']:.2f}%")
    print(f"Number of samples: {volume_stats['n_samples']}")


def save_analysis_to_file(analysis_results, volume_stats, filename='analysis_results.txt'):
    """
    Save analysis results to a text file.
    
    Parameters:
    -----------
    analysis_results : dict
        Results from calculate_water_density function
    volume_stats : dict
        Results from analyze_volume_fluctuations function
    filename : str
        Output filename
        
    Returns:
    --------
    str
        Path to the saved file
    """
    with open(filename, 'w', encoding='utf-8') as f:
        f.write("OpenMM Water Density Simulation Analysis\n")
        f.write("="*50 + "\n\n")
        
        f.write("Density Analysis:\n")
        f.write(f"Number of molecules: {analysis_results['n_molecules']}\n")
        f.write(f"Average volume: {analysis_results['avg_volume_nm3']:.6f} ± {analysis_results['std_volume_nm3']:.6f} nm³\n")
        f.write(f"Calculated density: {analysis_results['density']:.6f} g/cm³\n")
        f.write(f"Relative error: {calculate_density_error(analysis_results['density']):.3f}%\n\n")
        
        f.write("Volume Statistics:\n")
        f.write(f"Mean: {volume_stats['mean']:.6f} nm³\n")
        f.write(f"Std Dev: {volume_stats['std']:.6f} nm³\n")
        f.write(f"Min: {volume_stats['min']:.6f} nm³\n")
        f.write(f"Max: {volume_stats['max']:.6f} nm³\n")
        f.write(f"Range: {volume_stats['range']:.6f} nm³\n")
        f.write(f"CV: {volume_stats['cv']:.3f}%\n")
        f.write(f"Samples: {volume_stats['n_samples']}\n")
    
    print(f"Analysis results saved to '{filename}'")
    return filename


def complete_analysis(volumes, n_molecules, reference_density=1.0, 
                     save_to_file=True, output_file='analysis_results.txt'):
    """
    Perform complete analysis of simulation results.
    
    Parameters:
    -----------
    volumes : list
        List of volumes in nm³
    n_molecules : int
        Number of water molecules
    reference_density : float
        Reference density for comparison
    save_to_file : bool
        Whether to save results to file
    output_file : str
        Output filename
        
    Returns:
    --------
    dict
        Complete analysis results
    """
    # Calculate density
    density_results = calculate_water_density(volumes, n_molecules)
    
    # Analyze volume fluctuations
    volume_stats = analyze_volume_fluctuations(volumes)
    
    # Print results
    print_analysis_results(density_results, reference_density)
    print_volume_statistics(volume_stats)
    
    # Save to file if requested
    if save_to_file:
        save_analysis_to_file(density_results, volume_stats, output_file)
    
    # Calculate compressibility
    compressibility = calculate_compressibility(volumes)
    
    return {
        'density_results': density_results,
        'volume_stats': volume_stats,
        'compressibility': compressibility,
        'relative_error': calculate_density_error(density_results['density'], reference_density)
    }