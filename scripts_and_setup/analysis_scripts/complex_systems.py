from kimmdy_paper_theme import *
import matplotlib
from matplotlib.colors import Normalize, LinearSegmentedColormap, to_rgb
from matplotlib.cm import ScalarMappable
from scipy.stats import gaussian_kde
from matplotlib.patches import Patch
import os
import numpy as np

matplotlib.use("TkAgg")


def mimic_symmetric_alpha_colormap(hex_color, steps=256, lightness=0.5, darkness=0.5):
    base_rgb = np.array(to_rgb(hex_color))
    white = np.ones(3)
    black = np.zeros(3)

    light_rgb = white * (1 - lightness) + base_rgb * lightness
    dark_rgb = black * (1 - darkness) + base_rgb * darkness

    half_steps = steps // 2
    colors_light = np.linspace(light_rgb, base_rgb, half_steps, endpoint=False)
    colors_dark = np.linspace(base_rgb, dark_rgb, steps - half_steps)
    colors = np.vstack([colors_light, colors_dark])

    return LinearSegmentedColormap.from_list("symmetric_mimic_alpha_cmap", colors)


def kde_with_mass_levels(x_in, y_in, mass_levels_in, grid_x, grid_y):
    kde = gaussian_kde(np.vstack([x_in, y_in]))
    x_cur, y_cur = np.meshgrid(grid_x, grid_y)
    z_cur = kde(np.vstack([x_cur.ravel(), y_cur.ravel()])).reshape(x_cur.shape)

    z_sorted = np.sort(z_cur.ravel())
    cdf = np.cumsum(z_sorted)
    cdf /= cdf[-1]

    thresholds = []
    for m in mass_levels_in:
        idx = np.searchsorted(cdf, 1 - m)
        thresholds.append(z_sorted[idx])
    return x_cur, y_cur, z_cur, thresholds


np.random.seed(41)
plot_config = default_plot_config
apply_plot_config(plot_config=plot_config)
init_roboto_font()
plot_colors = plot_colors

green_cmap = mimic_symmetric_alpha_colormap(plot_colors["HITS_GREEN_LIGHT"], lightness=0.6, darkness=0.6)
magenta_cmap = mimic_symmetric_alpha_colormap(plot_colors["HITS_MAGENTA_LIGHT"], lightness=0.6, darkness=0.6)
yellow_cmap = mimic_symmetric_alpha_colormap(plot_colors["HITS_YELLOW_LIGHT"], lightness=0.6, darkness=0.6)
blue_cmap = mimic_symmetric_alpha_colormap(plot_colors["HITS_CYAN_LIGHT"], lightness=0.6, darkness=0.6)

data_dir = "../../data/data_KIMMDY/"
graphics_dir = "../../graphics/"

samples = set([a.split("R")[0].strip("_") for a in os.listdir(data_dir)])
quantum_yields = {sample: [] for sample in samples if not sample.startswith("ds_both")}
average_rates = {sample: [] for sample in samples if not sample.startswith("ds_both")}
distances = {sample: [] for sample in samples}
angles = {sample: [] for sample in samples}
quantum_yields["ds_both_post_left"] = []
quantum_yields["ds_both_post_right"] = []
quantum_yields["ds_both_pre_left"] = []
quantum_yields["ds_both_pre_right"] = []
average_rates["ds_both_post_left"] = []
average_rates["ds_both_post_right"] = []
average_rates["ds_both_pre_left"] = []
average_rates["ds_both_pre_right"] = []

window = 100

# Read data
for filename in os.listdir(data_dir):
    cur_file = open(f"{data_dir}{filename}")
    sample_name = filename.split("R")[0].strip("_")
    residues = (-1, -1)
    for line in cur_file:
        if line.startswith("Residues"):
            residues = line.strip("\n")
            residues = [int(a) for a in residues.split(" ")[1:]]
            residues = (residues[0], residues[1])
            current_key = residues
            continue
        if line.startswith("Rates"):
            rates = [float(a) for a in line.strip("\n").split("[")[1].strip("] ").split(",")]
            average_rate = sum(rates) / len(rates)
            quantum_yield = 100 * len([a for a in rates if a > 0.5]) / len(rates)
            if not filename.startswith("ds_both"):
                quantum_yields[sample_name].append(quantum_yield)
                average_rates[sample_name].append(average_rate)
            else:
                if residues == (13, 14):
                    quantum_yields[f"{sample_name}_left"].append(quantum_yield)
                    average_rates[f"{sample_name}_left"].append(average_rate)
                elif residues == (20, 21):
                    quantum_yields[f"{sample_name}_right"].append(quantum_yield)
                    average_rates[f"{sample_name}_right"].append(average_rate)

            print(
                f"{filename[0: len(filename) - 4]}: phi = {round(quantum_yield, 3)}, average rate = {round(average_rate, 3)}")
        if line.startswith("Distances"):
            distances_cur = [float(a) for a in line.strip("\n").split("[")[1].strip("] ").split(",")]
            distances[sample_name] += distances_cur
        if line.startswith("Angles"):
            angles_cur = [float(a) for a in line.strip("\n").split("[")[1].strip("] ").split(",")]
            angles[sample_name] += angles_cur

# Average both sites in crossovers
yields_cross = quantum_yields["crossover"]
quantum_yields["crossover"] = [(yields_cross[i] + yields_cross[i + 1]) / 2 for i in range(0, 12, 2)]
rates_cross = average_rates["crossover"]
average_rates["crossover"] = [(rates_cross[i] + rates_cross[i + 1]) / 2 for i in range(0, 12, 2)]

# Plotting quantum yields for origami variations
labels = ["ds_center", "nicked", "overhangs", "crossover"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05
quantum_yield_means = [np.mean(quantum_yields[l_i]) for l_i in labels]
quantum_yield_stds = [np.std(quantum_yields[l_i]) for l_i in labels]

group_colors = [plot_colors["HITS_MAGENTA_LIGHT"], plot_colors["HITS_YELLOW_LIGHT"], plot_colors["HITS_CYAN_LIGHT"],
                plot_colors["HITS_GREEN_LIGHT"]]
_, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], quantum_yield_means[i], width,
            yerr=quantum_yield_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Predicted Quantum Yield [%]")
ax1.set_ylim(bottom=0)
ax1.set_xticks(x)
labels_x = ["ds center", "nicked", "overhangs", "crossover"]
ax1.set_xticklabels(labels_x)
plt.savefig(f"{graphics_dir}ComparisonQuantumYieldOrigami.pdf")
plt.clf()

# Plotting rates origami
labels = ["ds_center", "nicked", "overhangs", "crossover"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05

average_rate_means = [np.mean(average_rates[l_i]) for l_i in labels]
average_rate_stds = [np.std(average_rates[l_i]) for l_i in labels]

group_colors = [plot_colors["HITS_MAGENTA_LIGHT"], plot_colors["HITS_YELLOW_LIGHT"], plot_colors["HITS_CYAN_LIGHT"],
                plot_colors["HITS_GREEN_LIGHT"]]
_, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], average_rate_means[i], width,
            yerr=average_rate_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Average reaction rate [ps$^{-1}$]")
ax1.set_ylim(bottom=0)
ax1.set_xticks(x)
labels_x = ["ds center", "nicked", "overhangs", "crossover"]
ax1.set_xticklabels(labels_x)
plt.savefig(f"{graphics_dir}ComparisonReactionRatesOrigami.pdf")
plt.clf()

# Plotting quantum yields consecutive
labels = ["ds_left", "ds_both_pre_left", "ds_both_post_left", "ds_right", "ds_both_pre_right", "ds_both_post_right"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05

quantum_yield_means = [np.mean(quantum_yields[l_i]) for l_i in labels]
quantum_yield_stds = [np.std(quantum_yields[l_i]) for l_i in labels]

group_colors = [plot_colors["MPG_grey_dark"] for i in range(0, 10)]
_, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], quantum_yield_means[i], width,
            yerr=quantum_yield_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Predicted Quantum Yield [%]")
ax1.set_xticks(x)
labels = ["left only", "left 1st", "left 2nd", "right only", "right 1st", "right 2nd"]
labels_x = labels
ax1.set_xticklabels(labels_x)
plt.savefig(f"{graphics_dir}ComparisonQuantumYieldConsecutive.pdf")
plt.clf()

# Plotting rates consecutive
labels = ["ds_left", "ds_both_pre_left", "ds_both_post_left", "ds_right", "ds_both_pre_right", "ds_both_post_right"]
x = np.arange(len(labels))
width = 0.3
spacing = 0.05

average_rate_means = [np.mean(average_rates[l_i]) for l_i in labels]
average_rate_stds = [np.std(average_rates[l_i]) for l_i in labels]

group_colors = [plot_colors["MPG_grey_dark"] for i in range(0, 10)]
_, ax1 = plt.subplots()

for i in range(len(labels)):
    ax1.bar(x[i], average_rate_means[i], width,
            yerr=average_rate_stds[i], color=group_colors[i], edgecolor="black",
            capsize=5)

ax1.set_ylabel("Average reaction rate [ps$^{-1}$]")
ax1.set_xticks(x)
labels = ["left only", "left 1st", "left 2nd", "right only", "right 1st", "right 2nd"]
labels_x = labels
ax1.set_xticklabels(labels_x)
plt.savefig(f"{graphics_dir}ComparisonReactionRatesConsecutive.pdf")
plt.clf()

# Random sample points for density plots
idx_cross = np.random.choice(len(distances["crossover"]), size=10000, replace=False)
idx_center = np.random.choice(len(distances["ds_center"]), size=10000, replace=False)
idx_nicked = np.random.choice(len(distances["nicked"]), size=10000, replace=False)
idx_overhangs = np.random.choice(len(distances["overhangs"]), size=10000, replace=False)

x_cross = np.array(distances["crossover"])[idx_cross]
y_cross = np.array(angles["crossover"])[idx_cross]
mask_cross = x_cross <= 5
x_cross, y_cross = x_cross[mask_cross], y_cross[mask_cross]

x_center = np.array(distances["ds_center"])[idx_center]
y_center = np.array(angles["ds_center"])[idx_center]
mask_tiny = x_center <= 5
x_center, y_center = x_center[mask_tiny], y_center[mask_tiny]

x_nicked = np.array(distances["nicked"])[idx_nicked]
y_nicked = np.array(angles["nicked"])[idx_nicked]
mask_nicked = x_nicked <= 5
x_nicked, y_nicked = x_nicked[mask_nicked], y_nicked[mask_nicked]

x_overhangs = np.array(distances["overhangs"])[idx_overhangs]
y_overhangs = np.array(angles["overhangs"])[idx_overhangs]
mask_overhangs = x_overhangs <= 5
x_overhangs, y_overhangs = x_overhangs[mask_overhangs], y_overhangs[mask_overhangs]

# KDE Comparison Plot

_, ax = plt.subplots()
mass_levels = np.linspace(0, 0.9, 10)
x_grid = np.linspace(0.2, 2.7, 1000)
y_grid = np.linspace(-205, 205, 1000)

X1, Y1, Z1, thresholds1 = kde_with_mass_levels(x_cross, y_cross, mass_levels, x_grid, y_grid)
levels1 = sorted(set(thresholds1))
Z1_masked = np.ma.masked_less(Z1, levels1[0])
cset1 = ax.contourf(X1, Y1, Z1_masked, levels=levels1, cmap=green_cmap)

X2, Y2, Z2, thresholds2 = kde_with_mass_levels(x_center, y_center, mass_levels, x_grid, y_grid)
levels2 = sorted(set(thresholds2))
Z2_masked = np.ma.masked_less(Z2, levels2[0])
cset2 = ax.contourf(X2, Y2, Z2_masked, levels=levels2, cmap=magenta_cmap, alpha=0.6)

cb_width = 0.016
cb_gap = 0.000

left1 = 0.91
left2 = left1 + cb_width + cb_gap

divider1 = plt.gcf().add_axes([left1, 0.15, cb_width, 0.7])
norm1 = Normalize(vmin=0, vmax=0.9)
sm1 = ScalarMappable(norm=norm1, cmap=green_cmap)
cb1 = plt.colorbar(sm1, cax=divider1, ticks=[])
cb1.ax.tick_params(size=0, width=0)
cb1.set_label("")

divider2 = plt.gcf().add_axes([left2, 0.15, cb_width, 0.7])
norm2 = Normalize(vmin=0, vmax=0.9)
sm2 = ScalarMappable(norm=norm2, cmap=magenta_cmap)
cb2 = plt.colorbar(sm2, cax=divider2, ticks=mass_levels)
cb2.set_ticklabels([f"{int(m * 100)}%" for m in mass_levels])
cb2.set_label("Normalized density")

ax.set_xlim(0.2, 2.7)
ax.set_ylim(-205, 205)
ax.set_xlabel("Distance [nm]")
ax.set_ylabel("Dihedral angle [deg]")
ax.set_facecolor("white")

legend_handles = [
    Patch(color=green_cmap(0.5), label="ds center"),
    Patch(color=magenta_cmap(0.5), label="crossover")
]

ax.legend(handles=legend_handles, loc="upper right", frameon=False)
plt.savefig(f"{graphics_dir}KDEComparison.pdf")
plt.clf()

# KDE Single Plots
# crossover
_, ax = plt.subplots()
mass_levels = np.linspace(0, 0.9, 10)
x_grid = np.linspace(0.2, 2.7, 1000)
y_grid = np.linspace(-205, 205, 1000)

X1, Y1, Z1, thresholds1 = kde_with_mass_levels(x_cross, y_cross, mass_levels, x_grid, y_grid)
levels1 = sorted(set(thresholds1))
Z1_masked = np.ma.masked_less(Z1, levels1[0])
_ = ax.contourf(X1, Y1, Z1_masked, levels=levels1, cmap=green_cmap)

cb_width = 0.016
left1 = 0.91
divider1 = plt.gcf().add_axes([left1, 0.15, cb_width, 0.7])
norm1 = Normalize(vmin=0, vmax=0.9)
sm1 = ScalarMappable(norm=norm1, cmap=green_cmap)
cb1 = plt.colorbar(sm1, cax=divider1, ticks=mass_levels)
cb1.set_ticklabels([f"{int(m * 100)}%" for m in mass_levels])
cb1.set_label("Normalized density")

ax.set_xlim(0.2, 2.7)
ax.set_ylim(-205, 205)
ax.set_xlabel("Distance [nm]")
ax.set_ylabel("Dihedral angle [deg]")
ax.set_facecolor("white")

plt.savefig(f"{graphics_dir}KDECrossover.pdf")
plt.clf()

# overhangs
_, ax = plt.subplots()
mass_levels = np.linspace(0, 0.9, 10)
x_grid = np.linspace(0.2, 2.7, 1000)
y_grid = np.linspace(-205, 205, 1000)

X1, Y1, Z1, thresholds1 = kde_with_mass_levels(x_overhangs, y_overhangs, mass_levels, x_grid, y_grid)
levels1 = sorted(set(thresholds1))
Z1_masked = np.ma.masked_less(Z1, levels1[0])
_ = ax.contourf(X1, Y1, Z1_masked, levels=levels1, cmap=blue_cmap)

cb_width = 0.016
left1 = 0.91
divider1 = plt.gcf().add_axes([left1, 0.15, cb_width, 0.7])
norm1 = Normalize(vmin=0, vmax=0.9)
sm1 = ScalarMappable(norm=norm1, cmap=blue_cmap)
cb1 = plt.colorbar(sm1, cax=divider1, ticks=mass_levels)
cb1.set_ticklabels([f"{int(m * 100)}%" for m in mass_levels])
cb1.set_label("Normalized density")

ax.set_xlim(0.2, 2.7)
ax.set_ylim(-205, 205)
ax.set_xlabel("Distance [nm]")
ax.set_ylabel("Dihedral angle [deg]")
ax.set_facecolor("white")

plt.savefig(f"{graphics_dir}KDEOverhangs.pdf")
plt.clf()

# nicked
_, ax = plt.subplots()
mass_levels = np.linspace(0, 0.9, 10)
x_grid = np.linspace(0.2, 2.7, 1000)
y_grid = np.linspace(-205, 205, 1000)

X1, Y1, Z1, thresholds1 = kde_with_mass_levels(x_nicked, y_nicked, mass_levels, x_grid, y_grid)
levels1 = sorted(set(thresholds1))
Z1_masked = np.ma.masked_less(Z1, levels1[0])
_ = ax.contourf(X1, Y1, Z1_masked, levels=levels1, cmap=yellow_cmap)

cb_width = 0.016
left1 = 0.91
divider1 = plt.gcf().add_axes([left1, 0.15, cb_width, 0.7])
norm1 = Normalize(vmin=0, vmax=0.9)
sm1 = ScalarMappable(norm=norm1, cmap=yellow_cmap)
cb1 = plt.colorbar(sm1, cax=divider1, ticks=mass_levels)
cb1.set_ticklabels([f"{int(m * 100)}%" for m in mass_levels])
cb1.set_label("Normalized density")

ax.set_xlim(0.2, 2.7)
ax.set_ylim(-205, 205)
ax.set_xlabel("Distance [nm]")
ax.set_ylabel("Dihedral angle [deg]")
ax.set_facecolor("white")

plt.savefig(f"{graphics_dir}KDENicked.pdf")
plt.clf()

# center
_, ax = plt.subplots()
mass_levels = np.linspace(0, 0.9, 10)
x_grid = np.linspace(0.2, 2.7, 1000)
y_grid = np.linspace(-205, 205, 1000)

X1, Y1, Z1, thresholds1 = kde_with_mass_levels(x_center, y_center, mass_levels, x_grid, y_grid)
levels1 = sorted(set(thresholds1))
Z1_masked = np.ma.masked_less(Z1, levels1[0])
_ = ax.contourf(X1, Y1, Z1_masked, levels=levels1, cmap=magenta_cmap)

cb_width = 0.016
left1 = 0.91
divider1 = plt.gcf().add_axes([left1, 0.15, cb_width, 0.7])
norm1 = Normalize(vmin=0, vmax=0.9)
sm1 = ScalarMappable(norm=norm1, cmap=magenta_cmap)
cb1 = plt.colorbar(sm1, cax=divider1, ticks=mass_levels)
cb1.set_ticklabels([f"{int(m * 100)}%" for m in mass_levels])
cb1.set_label("Normalized density")

ax.set_xlim(0.2, 2.7)
ax.set_ylim(-205, 205)
ax.set_xlabel("Distance [nm]")
ax.set_ylabel("Dihedral angle [deg]")
ax.set_facecolor("white")

plt.savefig(f"{graphics_dir}KDECenter.pdf")
plt.clf()
