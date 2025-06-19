from src.ljts.boks import Box
import os
import math
import numpy as np
from src.timeseries.autocorr import compute_autocorrelation, estimate_autocorr_time
from src.timeseries.block_avg_analysis import block_average, plot_block_averaging
import matplotlib.pyplot as plt
import time
def main():
    starting_time = time.time()
    length_x = 5
    length_y = 40
    length_z = 5
    density_liquid = 0.73
    density_vapor = 0.02

    box = Box(length_x, length_y, length_z, density_liquid, density_vapor)
    os.makedirs("./data", exist_ok=True)

    for path in [
        "./data/config_init.xyz",
        "./data/config_final.xyz",
        "./data/trajectory.xyz",
        "./data/sim_log.txt",
        "./data/results_log.txt"
    ]:
        if os.path.exists(path):
            os.remove(path)

    box.export_xyz("./data/config_init.xyz", "Initial configuration")
    initial_energy = box.compute_potential_energy()
    print(f"Initial potential energy: {initial_energy:.6f}")

    max_displacement = 1.0 / 8
    temperature = 0.8
    total_sweeps = 10000
    equilibration_sweeps = 5000
    log_interval = 10
    trajectory_interval = 200
    zeta = 1.00001

    energies = []
    exp_energies_1 = []
    exp_energies_2 = []

    surface_area_0 = box.get_surface_area()
    delta_area_1 = surface_area_0 * (zeta - 1)

    for sweep in range(1, total_sweeps + 1):
        box.mc_sweep(max_displacement, temperature)
        if (sweep % 100 == 0):
            print(f'Sweep: {sweep}/{total_sweeps} | {(sweep/total_sweeps)*100:.2f} %')
        if sweep > equilibration_sweeps and sweep % log_interval == 0:
            energy_undistorted = box.compute_potential_energy()
            energy_distorted_1 = box.compute_distorted_energy(sx=zeta**0.5, sy=1/zeta, sz=zeta**0.5)
            energy_distorted_2 = box.compute_distorted_energy(sx=1/zeta**0.5, sy=zeta, sz=1/zeta**0.5)

            exp_1 = math.exp(min(700, max(-700, -(energy_distorted_1 - energy_undistorted) / temperature)))
            exp_2 = math.exp(min(700, max(-700, -(energy_distorted_2 - energy_undistorted) / temperature)))

            energies.append(energy_undistorted)
            exp_energies_1.append(exp_1)
            exp_energies_2.append(exp_2)

        if sweep % trajectory_interval == 0:
            box.append_xyz_frame("./data/trajectory.xyz", f"Sweep {sweep}")

    compute_autocorrelation(
        energies, 
        lags=100, 
        plot=True, 
        title="Autocorrelation of Potential Energy",
        filename="./data/plots/autocorr_potential_energy.png"
    )

    compute_autocorrelation(
        exp_energies_1, 
        lags=100, 
        plot=True, 
        title="Autocorrelation of Exp(-ΔU/T) Scenario 1",
        filename="./data/plots/autocorr_exp_energy.png"
    )
    compute_autocorrelation(
        exp_energies_2, 
        lags=100, 
        plot=True, 
        title="Autocorrelation of Exp(-ΔU/T) Scenario 2",
        filename="./data/plots/autocorr_exp_energy_2.png"
    )

    block_size_energies = estimate_autocorr_time(energies)
    num_blocks_energies = len(energies) // block_size_energies

    mean_Epot, error_Epot, _ = block_average(energies, num_blocks_energies)

    block_size_exp1 = estimate_autocorr_time(exp_energies_1)
    block_size_exp2 = estimate_autocorr_time(exp_energies_2)
    num_blocks_exp1 = len(exp_energies_1) // block_size_exp1
    num_blocks_exp2 = len(exp_energies_2) // block_size_exp2
    mean_exp1, error_exp1, _ = block_average(exp_energies_1, num_blocks_exp1)
    mean_exp2, error_exp2, _ = block_average(exp_energies_2, num_blocks_exp2)
    gamma = -temperature * np.log(mean_exp1) / delta_area_1
    gamma_error = abs((temperature / (mean_exp1 * delta_area_1)) * error_exp1)

    plot_block_averaging(
        energies, 
        max_blocks=len(energies) // estimate_autocorr_time(energies), 
        title="Block Averaging of Potential Energy",
        filename="./data/plots/block_avg_potential_energy.png"
    )

    plot_block_averaging(
        exp_energies_1, 
        max_blocks=len(energies) // estimate_autocorr_time(exp_energies_1), 
        title="Block Averaging of Exp(-ΔU/T) Scenario 1",
        filename="./data/plots/block_avg_exp_energy.png"
    )
    plot_block_averaging(
        exp_energies_2, 
        max_blocks=len(energies) // estimate_autocorr_time(exp_energies_2), 
        title="Block Averaging of Exp(-ΔU/T) Scenario 1",
        filename="./data/plots/block_avg_exp_energy.png"
    )
    stop_time = time.time()
    run_time = (stop_time-starting_time)/60
    
    with open("./data/results_log.txt", "w", encoding="utf-8") as f:
        f.write(f"Initial potential energy: {initial_energy:.6f}\n")
        f.write(f"Optimal block size (energies): {block_size_energies}\n")
        f.write(f"<Epot> = {mean_Epot:.6f} ± {error_Epot:.6f}\n")
        f.write(f"<exp(-ΔU/T)> scenario 1 = {mean_exp1:.6f} ± {error_exp1:.6f}\n")
        f.write(f"<exp(-ΔU/T)> scenario 2 = {mean_exp2:.6f} ± {error_exp2:.6f}\n")
        f.write(f"Optimal block size (exp(-ΔU/T)): {block_size_exp1}\n")
        f.write(f"Surface tension γ = {gamma:.6f} ± {gamma_error:.6f}\n")
        f.write(f"Runtime (mins): {run_time:.2f}\n")
if __name__ == "__main__":
    main()
