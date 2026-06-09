import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
from scipy import special
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

# ============ Comparison of permeation as a function of time =================
# ========================= Diffusion limited =================================
num_summation_terms = 1000
#parameters listed/calculated from determine_limiter.py file
N_o = 6.11497e28 / 1000**3 #lattice density at/mm^3
lambdaa = 3.198e-10 * 1000 # lattice parameter (mm) || lambda is a python keyword
nu = 1e13 * 60 # Debye frequency (1/min)
rho = 0.01  # trapping site fraction
D_o = 7.44e-8 * 60 * 1000**2 # diffusivity pre-exponential (m^2/s)
Ed = 0.13  # diffusion activation energy (eV)
k = 8.617333262e-5  # Boltzmann's constant(eV/K)
T = 1500  # temperature (K)
epsilon = 1.2  # epsilon: trap energy (eV)
c = 0.0001  # dissolved gas atom fraction
zeta = ((lambdaa**2) * nu * np.exp((Ed - epsilon) / (k * T)) / (rho * D_o)) + (c / rho)
D = 2.72e-8 * 1000**2 * 60# diffusivity (mm^2/min)
D_eff = D / (1 + (1 / zeta))  # Effective diffusivity (mm^2/min)
l = 2.5e-5 * 1000  # slab thickness (mm)
c_o = c
tau_be = l**2 / (2 * (np.pi) ** 2 * D_eff)

def summation_term(num_terms, time):
    sum = 0.0 
    for m in range(1, num_terms):
        sum += (-1) ** m * np.exp(-1 * m**2 * time / (2 * tau_be))
    return sum

csv_folder = "./W-ver-1d-trapping_out.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_time = tmap_sol["time"]
tmap_perm = tmap_sol["scaled_outflux"]
idx = np.where(tmap_time >= 0.000001)[0][0]
c_o = c
analytical_time = tmap_time
Jp = (
     N_o * (c_o * D / l) * (1 + 2 * summation_term(num_summation_terms, analytical_time))
 )
# Plot figure for verification
x_max = 0.6
y_max = 1.2 * np.max(Jp)
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

analytical_permeation = Jp
ax.plot(
    analytical_time, analytical_permeation, label=r"Analytical", c="k", linestyle="--"
)
ax.plot(tmap_time, tmap_perm, label=r"TMAP8", c="tab:gray")

ax.set_xlabel("Time (min)")
ax.set_ylabel("Permeation (atom/mm$^2$min)")
ax.legend(loc="best")
ax.set_xlim(left=0, right=x_max)
ax.set_ylim(bottom=0, top = y_max)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
RMSE = np.sqrt(np.mean((tmap_perm - analytical_permeation)[idx:] ** 2))
RMSPE = RMSE * 100 / np.mean(analytical_permeation[idx:])
ax.text(0.65*x_max, 0.05*y_max, "RMSPE = %.2f " % RMSPE + "%", fontweight="bold")
ax.minorticks_on()
plt.savefig("tungsten-ver-1d_comparison_diffusion.png", bbox_inches="tight", dpi=300)
plt.close(fig)
# # ========================= Trapping limited ====================================
tau_bd = l**2 * rho / (2 * c_o * D)  # breakthrough time
analytical_time = [tau_bd, tau_bd]
steady_state_permeation = N_o * (c_o * D / l) 
analytical_sol = [0, steady_state_permeation]
csv_folder = "./W-ver-1d-trapping_out.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_time = tmap_sol["time"]
tmap_perm = tmap_sol["scaled_outflux"]
tmap_min_trapped = tmap_sol["min_trapped"]

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

y_max1 = 1.2 * steady_state_permeation
ax.plot(tmap_time, tmap_perm, label=r"TMAP8", c="tab:gray")
ax.plot(
    analytical_time,
    analytical_sol,
    label=r"Analytical breakthrough time",
    c="k",
    linestyle="--",
)
ax.axhline(
    y=steady_state_permeation,
    label=r"Analytical steady state",
    color="k",
    linestyle="-."
)
print(tmap_time[np.argmin(np.abs(tmap_perm - steady_state_permeation))])

ax.set_xlabel("Time (min)")
ax.set_ylabel("Permeation (atom/mm$^2$min)")
ax.legend(loc="best")
ax.set_xlim(left=0, right=x_max)
ax.set_ylim(bottom=0, top = y_max1)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig("tungsten-ver-1d_comparison_trapping.png", bbox_inches="tight", dpi=300)
plt.close(fig)
