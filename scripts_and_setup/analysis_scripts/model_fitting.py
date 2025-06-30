from kimmdy_paper_theme import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math

matplotlib.use("TkAgg")


def calculate_rate(k_1_in, k_2_in, distance, angle):
    optimal_distance = 0.157177  # nm
    optimal_angle = 16.743651884789273  # deg
    return math.exp(-(k_1_in * abs(distance - optimal_distance) + k_2_in * abs(angle - optimal_angle)))


plot_config = default_plot_config
apply_plot_config(plot_config=plot_config)
init_roboto_font()
plot_colors = plot_colors

grid_k1 = np.linspace(0, 5, 1000)
grid_k2 = np.linspace(0, 1, 1000)

dT20_quantum_yield_diffs_perT = np.load("../../data/data_dT20/dT20_1000x1000_perT.npy")
plt.figure()
plt.imshow(dT20_quantum_yield_diffs_perT, extent=[grid_k2.min(), grid_k2.max(), grid_k1.min(), grid_k1.max()],
           origin='lower', aspect='auto', cmap='viridis')
cb = plt.colorbar()
cb.set_label(label=r'|$\phi-\phi_{\mathrm{exp}}|$')
plt.xlabel(r'$k_2$')
plt.ylabel(r'$k_1$')
# plt.title(r'dT$_{20}$, normlization per T')
plt.tight_layout()
plt.savefig("../../graphics/grid_dT20_T.pdf")
plt.clf()

TdT_quantum_yield_diffs = np.load("../../data/data_TdT/TdT_1000x1000.npy")
plt.figure()
plt.imshow(TdT_quantum_yield_diffs, extent=[grid_k2.min(), grid_k2.max(), grid_k1.min(), grid_k1.max()],
           origin='lower', aspect='auto', cmap='viridis')
cb = plt.colorbar()
cb.set_label(label=r'|$\phi-\phi_{\mathrm{exp}}|$')
plt.xlabel(r'$k_2$')
plt.ylabel(r'$k_1$')
# plt.title('TdT')
plt.tight_layout()
plt.savefig("../../graphics/grid_TdT.pdf")
plt.clf()

max_TdT_dT20_perT = np.maximum(dT20_quantum_yield_diffs_perT, TdT_quantum_yield_diffs)
min_idx = np.unravel_index(np.argmin(max_TdT_dT20_perT), max_TdT_dT20_perT.shape)

print("Minimum value for T fit is at:", np.min(max_TdT_dT20_perT), grid_k1[min_idx[0]], grid_k2[min_idx[1]])
print("Minimum values in matrix", dT20_quantum_yield_diffs_perT[min_idx[0], [min_idx[1]]],
      TdT_quantum_yield_diffs[min_idx[0], [min_idx[1]]])

plt.figure()
plt.imshow(max_TdT_dT20_perT, extent=[grid_k2.min(), grid_k2.max(), grid_k1.min(), grid_k1.max()],
           origin='lower', aspect='auto', cmap='viridis')
cb = plt.colorbar()
cb.set_label(label=r'$\mathrm{max}(|\Delta \phi_{\mathrm{TdT}}, \Delta \phi_{\mathrm{dT}_{20}})$')
plt.xlabel(r'$k_2$')
plt.ylabel(r'$k_1$')
# plt.title(r"Max $\Delta \phi$ TdT and dT$_{20}$, normlization per T")
plt.plot(grid_k2[min_idx[1]], grid_k1[min_idx[0]], 'rx')
plt.tight_layout()
plt.savefig("../../graphics/grid_max_T.pdf")
plt.clf()

print("Law:", calculate_rate(grid_k1[min_idx[0]], grid_k2[min_idx[1]], 0.37, 48.2))
print("Johnson:", calculate_rate(grid_k1[min_idx[0]], grid_k2[min_idx[1]], 0.34, 27.0))
print("McCullagh:", calculate_rate(grid_k1[min_idx[0]], grid_k2[min_idx[1]], 0.352, 40))
