"""
Run the water density workflow across a temperature ramp.

This script reuses the single-temperature routine defined in
`water_density_at_temperature.py` and sweeps from 280 K up to
310 K in 5 K increments, printing a concise summary after each
simulation and a consolidated table at the end.
"""

import matplotlib.pyplot as plt
from openmm import unit

from water_density_at_temperature import density_at_temperature


def run_temperature_ramp(
    start_k: int = 260,
    stop_k: int = 310,
    step_k: int = 2,
) -> list[dict]:
    """
    Execute the density workflow for a sequence of temperatures.

    Returns a list of dictionaries containing the temperature (K)
    and the resulting density (g/cm^3) for convenience.
    """
    results = []
    for temp_k in range(start_k, stop_k + 1, step_k):
        temperature = temp_k * unit.kelvin
        print("\n" + "=" * 60)
        print(f"Running water density simulation at {temp_k} K")
        print("=" * 60)

        analysis = density_at_temperature(temperature=temperature)
        density_value = analysis['density_results']['density']
        rel_error = analysis['relative_error']

        results.append(
            {
                'temperature_K': temp_k,
                'density_g_per_cm3': density_value,
                'relative_error_percent': rel_error,
            }
        )

        print(
            f"Completed {temp_k} K "
            f"→ density {density_value:.4f} g/cm³ "
            f"(error {rel_error:.2f}%)"
        )

    return results


def plot_density_results(
    ramp_results: list[dict],
    output_file: str = "density_vs_temperature.png"
):
    """
    Plot density vs. temperature for the ramp results and save to disk.
    """
    if not ramp_results:
        print("No ramp results available to plot.")
        return None

    if plt is None:
        print("matplotlib is not installed; skipping density plot.")
        return None

    temperatures = [entry['temperature_K'] for entry in ramp_results]
    densities = [entry['density_g_per_cm3'] for entry in ramp_results]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(temperatures, densities, marker='o', linestyle='-', color='tab:blue')
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Density (g/cm³)")
    ax.set_title("Water Density vs Temperature")
    ax.grid(True, linestyle='--', alpha=0.5)

    fig.tight_layout()
    fig.savefig(output_file, dpi=300)
    plt.close(fig)

    print(f"Density plot saved to '{output_file}'")
    return output_file


def save_ramp_summary(
    ramp_results: list[dict],
    filename: str = "temperature_ramp_summary.txt"
):
    """
    Persist the ramp summary table to a text file.
    """
    if not ramp_results:
        print("No ramp results available to save.")
        return None

    lines = [
        "Temperature Ramp Summary",
        "#" * 60,
    ]
    for item in ramp_results:
        lines.append(
            f"{item['temperature_K']:>3} K | "
            f"density {item['density_g_per_cm3']:.4f} g/cm³ | "
            f"error {item['relative_error_percent']:.2f}%"
        )

    with open(filename, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))

    print(f"Ramp summary written to '{filename}'")
    return filename


def main():
    ramp_results = run_temperature_ramp()

    print("\n" + "#" * 60)
    print("Temperature Ramp Summary")
    print("#" * 60)
    for item in ramp_results:
        print(
            f"{item['temperature_K']:>3} K | "
            f"density {item['density_g_per_cm3']:.4f} g/cm³ | "
            f"error {item['relative_error_percent']:.2f}%"
        )

    plot_density_results(ramp_results)
    save_ramp_summary(ramp_results)


if __name__ == "__main__":
    main()
