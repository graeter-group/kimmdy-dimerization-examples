from kimmdy_paper_theme import *
import matplotlib
import math
import numpy as np

matplotlib.use("TkAgg")


def read_angle_file(file_name_in):
    file = open(file_name_in)
    data = False
    data_angles = {}

    for line_cur in file:
        if "@TYPE xy" in line_cur and not data:
            data = True
            continue
        if data:
            line_cur = [a for a in line_cur.strip("\n").split(" ") if a]
            time = int(float(line_cur[0]))  # in ps
            values = [float(a) for a in line_cur[1:]]
            data_angles[time] = values
    return data_angles


def read_distance_file(file_name_in):
    file = open(file_name_in)
    data = False
    data_distances = {}

    for line_cur in file:
        if "@TYPE xy" in line_cur and not data:
            data = True
            continue
        if data:
            line_cur = [a for a in line_cur.strip("\n").split(" ") if a]
            time = int(float(line_cur[0]))  # in ps
            residues = int((len(line_cur) - 1) / 2 + 1)
            values = [(float(line_cur[2 * j - 1]) + float(line_cur[2 * j])) / 2 for j in range(1, residues)]  # in nm
            data_distances[time] = values
    return data_distances


def calculate_rate(distance_in, angle_in):
    optimal_distance = 0.157177  # nm
    optimal_angle = 16.743651884789273  # deg
    k_1 = 2.017017017017017
    k_2 = 0.03003003003003003
    return math.exp(-(k_1 * abs(distance_in - optimal_distance) + k_2 * abs(angle_in - optimal_angle)))


plot_config = default_plot_config
apply_plot_config(plot_config=plot_config)
init_roboto_font()
plot_colors = plot_colors

data_dir = "../../data/data_TdT_KIMMDY/"
graphics_dir = "../../graphics/"

runs = [2, 3]

window = 100

# Read Dihedrals
for run in runs:
    for base in range(1, 2):
        cur_file = open(f"{data_dir}TdT_R{run}_Dihedral{base}.xvg")
        times = []
        dihedrals = []
        for line in cur_file:
            line = [float(a) for a in line.strip("\n").split(" ") if a != ""]
            times.append(line[0] / 1000)
            dihedrals.append(line[1])
        running_average = []
        for idx in range(0 + int(window / 2), len(dihedrals) - int(window / 2)):
            running_average.append(sum(dihedrals[int(idx - window / 2):int(idx + window / 2)]) / window)

        if run == 2:
            plt.plot(times[0 + int(window / 2): len(dihedrals) - int(window / 2)], running_average,
                     label=f"Simulation 1", color=plot_colors["HITS_MAGENTA_LIGHT"])
        else:
            plt.plot(times[0 + int(window / 2): len(dihedrals) - int(window / 2)], running_average,
                     label=f"Simulation 2", color=plot_colors["HITS_GREEN_LIGHT"])

plt.ylabel(r"$\chi$ (deg)")
plt.xlabel(r"Simulation time [ns]")

plt.legend()
plt.tight_layout()
plt.savefig(f"{graphics_dir}DihedralComparisonTdT.pdf")
plt.clf()

quantum_yields = {"Anti-Anti": [], "Anti-Syn": [], "dT20": []}
average_rates = {"Anti-Anti": [], "Anti-Syn": [], "dT20": []}

runs = [1, 2, 3]
flip_frames = [23265, 4590, None]

# Read data from KIMMDY Runs [Run 1, Run 2 are Anti-Syn, Run 3 is Anti-Anti]
for run, flip_frame in zip(runs, flip_frames):
    cur_file = open(f"{data_dir}TdT_R{run}_reaction_rates.csv")

    for line in cur_file:
        if not line.startswith("Rates"):
            continue
        else:
            rates = [float(a) for a in line.strip("\n").split("[")[1].strip("] ").split(",")]
            if flip_frame:
                rates = rates[flip_frame:]
            average_rate = sum(rates) / len(rates)
            quantum_yield = 100 * len([a for a in rates if a > 0.5]) / len(rates)
            plt.xlabel(r"Rate [ps$^{-1}$]")
            plt.ylabel("# of Frames")
            if run != 3:
                quantum_yields["Anti-Syn"].append(quantum_yield)
                average_rates["Anti-Syn"].append(average_rate)
                plt.hist(rates, bins=50, color=plot_colors["HITS_MAGENTA_LIGHT"])
            else:
                quantum_yields["Anti-Anti"].append(quantum_yield)
                average_rates["Anti-Anti"].append(average_rate)
                plt.hist(rates, bins=50, color=plot_colors["HITS_GREEN_LIGHT"])

            plt.tight_layout()

            plt.savefig(f"{graphics_dir}HistogramRatesTdT_R{run}.pdf")
            plt.clf()

data_dir = "../../data/data_TdT/"

# Read data from non KIMMDY Runs TdT
old_runs = [1, 2, 3]
for run in old_runs:
    angles = list(read_angle_file(f"{data_dir}TdT_R{run}_MD_Angles.xvg").values())
    distances = list(read_distance_file(f"{data_dir}TdT_R{run}_MD_Distances.xvg").values())
    rates = []
    for angle, distance in zip(angles, distances):
        rates.append(calculate_rate(distance[0], angle[0]))
    average_rate = sum(rates) / len(rates)
    quantum_yield = 100 * len([a for a in rates if a > 0.5]) / len(rates)
    quantum_yields["Anti-Anti"].append(quantum_yield)
    average_rates["Anti-Anti"].append(average_rate)

data_dir = "../../data/data_dT20/"

# Read data from non KIMMDY Runs dT20
old_runs = [1, 2, 3]
for run in old_runs:
    angles = list(read_angle_file(f"{data_dir}dT20_R{run}_MD_Angles.xvg").values())
    distances = list(read_distance_file(f"{data_dir}dT20_R{run}_MD_Distances.xvg").values())
    rates = []
    for angle, distance in zip(angles, distances):
        for angle_sing, distance_sing in zip(angle, distance):
            rates.append(calculate_rate(distance_sing, angle_sing))
    average_rate = sum(rates) / len(rates)
    quantum_yield = 2 * 100 * len([a for a in rates if a > 0.5]) / len(rates)
    quantum_yields["dT20"].append(quantum_yield)
    average_rates["dT20"].append(average_rate)

# Quantum yields
labels = ["Anti-Anti", "Anti-Syn", "dT20"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05

quantum_yield_means = [np.mean(quantum_yields[l_i]) for l_i in labels]
quantum_yield_stds = [np.std(quantum_yields[l_i]) for l_i in labels]
average_rate_means = [np.mean(average_rates[l_i]) for l_i in labels]
average_rate_stds = [np.std(average_rates[l_i]) for l_i in labels]

group_colors = [plot_colors["HITS_GREEN_LIGHT"], plot_colors["HITS_MAGENTA_LIGHT"], plot_colors["HITS_CYAN_LIGHT"]]

fig, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], quantum_yield_means[i], width,
            yerr=quantum_yield_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Predicted Quantum Yield [%]")
ax1.set_xticks(x)
labels_x = [r"TdT $\mathit{anti}$-$\mathit{anti}$", r"TdT $\mathit{anti}$-$\mathit{syn}$", "dT20"]
ax1.set_xticklabels(labels_x)
fig.tight_layout()
plt.savefig(f"{graphics_dir}ComparisonQuantumYieldTdT.pdf")

# Average Rates
labels = ["Anti-Anti", "Anti-Syn", "dT20"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05

average_rate_means = [np.mean(average_rates[l_i]) for l_i in labels]
average_rate_stds = [np.std(average_rates[l_i]) for l_i in labels]

group_colors = [plot_colors["HITS_GREEN_LIGHT"], plot_colors["HITS_MAGENTA_LIGHT"], plot_colors["HITS_CYAN_LIGHT"]]

fig, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], average_rate_means[i], width,
            yerr=average_rate_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Average reaction rate [ps$^{-1}$]")
ax1.set_xticks(x)
labels_x = [r"TdT $\mathit{anti}$-$\mathit{anti}$", r"TdT $\mathit{anti}$-$\mathit{syn}$", "dT20"]
ax1.set_xticklabels(labels_x)
fig.tight_layout()
plt.savefig(f"{graphics_dir}ComparisonReactionRatesTdT.pdf")
