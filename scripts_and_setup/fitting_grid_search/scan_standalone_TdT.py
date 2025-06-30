from read_data import *
import math
import sys
import numpy as np


def calculate_rate(k_1_in, k_2_in, distance, angle):
    optimal_distance = 0.157177  # nm
    optimal_angle = 16.743651884789273  # deg
    return math.exp(-(k_1_in * abs(distance - optimal_distance) + k_2_in * abs(angle - optimal_angle)))


runs = 3
quantum_yield_exp = 0.013
n_T = 2

# Read distances
distances_dT20 = {}
for i in range(1, runs+1):
    distances_dT20[f"Run {i}"] = read_distance_file(f"../../data/data_TdT/TdT_R{i}_MD_Distances.xvg")

# Read dihedrals
dihedrals_dT20 = {}
for i in range(1, runs+1):
    dihedrals_dT20[f"Run {i}"] = read_angle_file(f"../../data/data_TdT/TdT_R{i}_MD_Angles.xvg")

start_k1 = float(sys.argv[1])
end_k1 = float(sys.argv[2])
spacing_k1 = int(sys.argv[3])
start_k2 = float(sys.argv[4])
end_k2 = float(sys.argv[5])
spacing_k2 = int(sys.argv[6])
job_index = sys.argv[7]

grid_k1 = np.linspace(start_k1, end_k1, spacing_k1)
grid_k2 = np.linspace(start_k2, end_k2, spacing_k2)
quantum_yield_diffs = np.zeros((len(grid_k1), len(grid_k2)))

for i, k_1 in enumerate(grid_k1):
    for j, k_2 in enumerate(grid_k2):
        successful = 0
        times = []
        for run in range(1, runs+1):
            distances_cur = distances_dT20[f"Run {run}"]
            dihedrals_cur = dihedrals_dT20[f"Run {run}"]
            times = sorted(list(distances_cur.keys()))
            for time in distances_cur.keys():
                distances_cur_time = distances_cur[time]
                dihedrals_cur_time = dihedrals_cur[time]
                for distance_cur, dihedral_cur in zip(distances_cur_time, dihedrals_cur_time):
                    rate_cur = calculate_rate(k_1, k_2, distance_cur, abs(dihedral_cur))
                    if rate_cur > 0.5:
                        successful += 1
        quantum_yield = min(quantum_yield_exp, abs(successful*(1/(n_T-1))*1/(runs*len(times))-quantum_yield_exp))
        quantum_yield_diffs[i, j] = quantum_yield

np.save(f"{job_index}.npy", quantum_yield_diffs)



