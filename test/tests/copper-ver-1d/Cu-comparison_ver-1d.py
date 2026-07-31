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
F = 96485.332 #C/mol, Faraday's constant
T = 1500  # temperature (K)
R = 8.31446261815324 #J/molK, gas constant
k = 8.617333262e-5  # Boltzmann's constant(eV/K)

lambdaa = 1e-4*(3.602472 + 6.99e-6*T + 1.45e-7*T**2)  # lattice parameter (mum) || lambda is a python keyword
N_o = 4/(lambdaa**3) #lattice density at/mum^3
nu = 1e13 # Debye frequency (1/s)
rho = 0.01  # trapping site fraction
D_o = 1.74e-6 * 1e12# diffusivity pre-exponential (mum^2/s)
Ed = 42000/F  # diffusion activation energy (eV)
Ed_J_mol = 42000 #diffusion activation energy (J/mol)
epsilon = 0.26  # epsilon: trap energy (eV)
c = 0.0001  # dissolved gas atom fraction
zeta = ((lambdaa**2) * nu * np.exp((Ed - epsilon) / (k * T)) / (rho * D_o)) + (c / rho)
D = D_o*np.exp(-Ed_J_mol/(R*T))# diffusivity (mum^2/s)
D_eff = D / (1 + (1 / zeta))  # Effective diffusivity (mum^2/min)
l = 50  # slab thickness (mum)
c_o = c
tau_be = l**2 / (2 * (np.pi) ** 2 * D_eff)

def summation_term(num_terms, time):
    sum = 0.0 
    for m in range(1, num_terms):
        sum += (-1) ** m * np.exp(-1 * m**2 * time / (2 * tau_be))
    return sum

csv_folder = "./Cu-ver-1d-trapping_out.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_time = tmap_sol["time"]
tmap_perm = tmap_sol["scaled_outflux"]
idx = np.where(tmap_time >= 1e-5)[0][0]

c_o = c
analytical_time = tmap_time
Jp = (
     N_o * (c_o * D / l) * (1 + 2 * summation_term(num_summation_terms, analytical_time))
 )
# Plot figure for verification
x_max = 0.2
y_max = 1.2 * np.max(Jp)
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

analytical_permeation = Jp
ax.plot(
    analytical_time, analytical_permeation, label=r"Analytical", c="k", linestyle="--"
)
ax.plot(tmap_time, tmap_perm, label=r"TMAP8", c="tab:gray")

ax.set_xlabel("Time (s)")
ax.set_ylabel("Permeation (atom/m$^2$s)")
ax.legend(loc="best")
ax.set_xlim(left=0, right=x_max)
ax.set_ylim(bottom=0, top = y_max)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
RMSE = np.sqrt(np.mean((tmap_perm - analytical_permeation)[idx:] ** 2))
RMSPE = RMSE * 100 / np.mean(analytical_permeation[idx:])
ax.text(0.65*x_max, 0.05*y_max, "RMSPE = %.2f " % RMSPE + "%", fontweight="bold")
ax.minorticks_on()
plt.savefig("copper-ver-1d_comparison_diffusion.png", bbox_inches="tight", dpi=300)
plt.close(fig)
# # ========================= Trapping limited ====================================
tau_bd = l**2 * rho / (2 * c_o * D)  # breakthrough time
analytical_time = [tau_bd, tau_bd]
steady_state_permeation = N_o * (c_o * D / l) 
print(steady_state_permeation)
analytical_sol = [0, steady_state_permeation]
csv_folder = "./Cu-ver-1d-trapping_out.csv"
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
    color="g",
    linestyle="-."
)
print(tmap_time[np.argmin(np.abs(tmap_perm - steady_state_permeation))])

ax.set_xlabel("Time (s)")
ax.set_ylabel("Permeation (atom/m$^2$s)")
ax.legend(loc="best")
ax.set_xlim(left=0, right=x_max)
ax.set_ylim(bottom=0, top = y_max1)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig("copper-ver-1d_comparison_trapping.png", bbox_inches="tight", dpi=300)
plt.close(fig)
