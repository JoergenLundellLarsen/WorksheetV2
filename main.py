from src.ljts.boks import Box
import matplotlib.pyplot as plt
import os
import math


def main():
    # Box dimensions and densities
    length_x = 5
    length_y = 40
    length_z = 5
    density_liquid = 0.73
    density_vapor = 0.02

    # Create simulation box
    box = Box(length_x, length_y, length_z, density_liquid, density_vapor)
    os.makedirs("./data", exist_ok=True)

    # Exportconfiguration
    box.export_xyz("./data/config_init.xyz", "Initial configuration")

    initial_energy = box.compute_potential_energy()
    print(f"Initial potential energy: {initial_energy:.6f}")

    num_molecules = box.get_molecule_count()
    max_displacement = 1.0 / 8
    temperature = 0.8
    total_sweeps = 1000
    equilibration_sweeps = 500
    log_interval = 100
    trajectory_interval = 200
    zeta = 1.00001

    sum_expansion_1 = 0.0
    sum_expansion_2 = 0.0
    sum_energy = 0.0
    sample_count = 0

    log_lines = []
    log_filename = "./data/sim_log.txt"

    # Original surface area
    surface_area_0 = box.get_surface_area()
    # Surface area differences due to distortions
    delta_area_1 = surface_area_0 * (zeta - 1)
    delta_area_2 = surface_area_0 * (1/zeta - 1)

    for sweep in range(1, total_sweeps + 1):
        box.mc_sweep(max_displacement, temperature)

        if sweep > equilibration_sweeps:
            energy_undistorted = box.compute_potential_energy()
            energy_distorted_1 = box.compute_distorted_energy(sx=zeta, sy=1/zeta, sz=zeta)
            energy_distorted_2 = box.compute_distorted_energy(sx=1/zeta, sy=zeta, sz=1/zeta)

            exponent_1 = max(-700, min(700, -(energy_distorted_1 - energy_undistorted) / temperature))
            exponent_2 = max(-700, min(700, -(energy_distorted_2 - energy_undistorted) / temperature))

            exp_1 = math.exp(exponent_1)
            exp_2 = math.exp(exponent_2)

            sum_expansion_1 += exp_1
            sum_expansion_2 += exp_2
            sum_energy += energy_undistorted
            sample_count += 1

        if sweep % trajectory_interval == 0:
            box.append_xyz_frame("./data/trajectory.xyz", f"Sweep {sweep}")

        if sweep % log_interval == 0 and sweep > equilibration_sweeps:
            avg_exp_1 = sum_expansion_1 / sample_count
            avg_exp_2 = sum_expansion_2 / sample_count
            avg_energy = sum_energy / sample_count

            gamma_1 = -temperature * math.log(avg_exp_1) / delta_area_1 if avg_exp_1 > 0 else float('inf')
            gamma_2 = -temperature * math.log(avg_exp_2) / delta_area_2 if avg_exp_2 > 0 else float('inf')

            log = (
                f"Sweep {sweep}\n"
                f" Current Epot: {energy_undistorted:.6f} | Avg: {avg_energy:.6f}\n"
                f" <exp(-ΔU/T)> s1: {avg_exp_1:.6f} → γ₁: {gamma_1:.6f}\n"
                f" <exp(-ΔU/T)> s2: {avg_exp_2:.6f} → γ₂: {gamma_2:.6f}\n"
                f"{'-'*40}"
            )
            print(log)
            log_lines.append(log)

    # Export final configuration and save logs
    box.export_xyz("./data/config_final.xyz", "Final configuration")
    with open(log_filename, "w", encoding="utf-8") as f:
        f.write("\n".join(log_lines))

    # Final measurement 
    final_energy = box.compute_potential_energy()
    scale_x, scale_y, scale_z = 1.1, 1.0, 1.1
    distorted_energy = box.compute_distorted_energy(scale_x, scale_y, scale_z)
    delta_energy = distorted_energy - final_energy
    print(f"\nΔU (energy change after distortion): {delta_energy:.6f}")

    surface_area_new = surface_area_0 * scale_x * scale_z
    delta_surface_area = surface_area_new - surface_area_0
    print(f"ΔA (surface area change): {delta_surface_area:.6f}")


if __name__ == "__main__":
    main()