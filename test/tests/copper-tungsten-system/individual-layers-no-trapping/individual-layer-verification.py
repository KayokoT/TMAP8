from ast import main
import string

from importlib_metadata import files
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
from scipy import special
from numpy import pi, sin, cos, tan, sqrt, exp
import os
import glob
import re
from collections import defaultdict
# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

"""
 This file performs verification checks on each run of Cu-W-Ar.i by plotting the conservation of mass over time. 
 It also compares the diffusion front to the analytical solution for Cu and W.
 It also creates plots the concentration over time and the concentration over distance for the Cu and W individual
 cases.
"""

# === Input Parameters ======
R = 8.31446261815324   #J/molK, gas constant
kb = 8.617333262e-5    #eV/K Boltzmann's constant
Cu_thickness = 50
W_thickness = 25
T = 900
Cu_D0 = 1.74e-6                    #m^2/s, diffusivity pre-exponential
Cu_Ed = 42000                                     #J/mol, activation energy
Cu_D = Cu_D0*exp(-Cu_Ed/(R*T)) * 1e12 #mum^2/s diffusivity
W_D0 = 7.44e-8                  #m^2/s, diffusivity pre-exponential [1]
W_Ed = 0.13                                      #eV, diffusion activation energy [1]
W_D = W_D0*exp(-W_Ed/(kb*T)) * 1e12 #mum^2/s, diffusivity

# =================================== Function ===============================
def read_csv_from_TMAP8(file_name, parameter_names):
    """Read simulation output columns from a CSV file into a numpy array
    From AUNF Mini-Canister Example Case

    Args:
        file_name (str): name of the CSV file in the gold directory
        parameter_names (list of str): column names to extract, in desired order

    Returns:
        ndarray: 2D array of shape (len(parameter_names), n_timesteps)
    """
    csv_folder = f"./{file_name}"
    simulation_data = pd.read_csv(csv_folder)
    return np.array([simulation_data[name] for name in parameter_names])

def compute_rmspe(simulated, reference):
    """Compute the Root Mean Square Percentage Error between two arrays

    Args:
        simulated (float, ndarray): simulated values
        reference (float, ndarray): reference values used as the denominator

    Returns:
        float: RMSPE in percent
    """
    RMSE = np.sqrt(np.mean((simulated - reference) ** 2))
    return RMSE * 100 / np.mean(reference)

def annotate_rmspe(simulated, reference, x_pos, y_pos):
    """Compute RMSPE and annotate it as bold text on the current matplotlib axes

    Args:
        simulated (float, ndarray): simulated values
        reference (float, ndarray): reference values used as the denominator
        x_pos (float): x-coordinate of the annotation
        y_pos (float): y-coordinate of the annotation
    """
    RMSPE = compute_rmspe(simulated, reference)
    plt.text(x_pos, y_pos, "RMSPE = %.2f %%" % RMSPE, fontweight="bold")

def plot_conservation_of_mass(t, flux, mass, t0, filename):
    """Plot the absolute percent difference between accumulated boundary flux
    and total mass as a single curve, verifying mass conservation. Data prior
    to t = t0 s is excluded so the metric is not dominated by early-time
    relative noise from the time-integrated flux postprocessor — analogous to
    the early-time slicing convention used for RMSPE in the ver-* cases.

    Args:
        t (float, ndarray): time array in seconds
        t0 (float) : time before which is excluded from metric
        flux (float, ndarray): accumulated boundary flux in mol H
        mass (float, ndarray): total H mass in domain in mol H
        filename (str): output PNG filename
    """
    t = np.asarray(t)
    flux = np.asarray(flux, dtype=float)
    mass = np.asarray(mass, dtype=float)
    idx = np.where(t > t0)[0][0]
    t = t[idx:]
    percent_diff = np.abs(flux[idx:] / mass[idx:] - 1.0) * 100.0

    fig = plt.figure(figsize=[6.5, 5.5])
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(t, percent_diff, c="tab:red")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Percent difference (%)")
    ax.set_xlim(left=t[0], right=t[-1])
    ax.set_ylim(bottom=0)
    plt.grid(which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.set_title("Conservation of Mass Check")
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig)

# ==================================== Copper only ================================================
# ================ Data Extraction ================
Cu_time, conc_Cu, Cu_exact_diff_length, Cu_sim_diff_length, Cu_mass_from_flux, Cu_mass_from_concentration = read_csv_from_TMAP8("copper-only_main.csv", 
                                                                                                                  ["time", "concentration_at_x_Cu",
                                                                                                                   "exact_diffusion_length", "simulated_diffusion_length",
                                                                                                                   "mass_from_flux", "mass_from_concentration"])
# ================ Plot Concentration over Time ================
fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(Cu_time, conc_Cu, color="mediumslateblue")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Concentration (at/μm³)")
ax.set_xlim(0, np.max(Cu_time))
ax.set_ylim(0, 1.2* np.max(conc_Cu))
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Copper-only concentration vs. time")
fig.savefig("copper-only-conc-vs-time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Plot Concentration over Distance ================
Cu_x, conc_Cu = read_csv_from_TMAP8(sorted(glob.glob("copper-only_vector_Cu_*.csv"))[-1], ["x", "C_M_Cu"])

fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(Cu_x, conc_Cu, color="darkred")
ax.set_xlabel("Distance (μm)")
ax.set_ylabel("Concentration (at/μm³)")
ax.set_xlim(0, np.max(Cu_x))
ax.set_ylim(0, 1.2*np.max(conc_Cu))
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Copper-only concentration vs. distance")
fig.savefig("copper-only-conc-vs-distance.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Diffusion Front Check ================
# Constant-pressure run: used only for the diffusion-front verification, whose
# analytical solution ell(t) = sqrt(pi*D*t) assumes a step Dirichlet BC.
fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(Cu_time, Cu_exact_diff_length, label = "Analytical", color="darkviolet")
ax.plot(Cu_time, Cu_sim_diff_length, label = "TMAP8", color="darkred")
annotate_rmspe(
    Cu_sim_diff_length,
    Cu_exact_diff_length,
    Cu_time[-1] / 2,
    Cu_exact_diff_length[-1] / 4,
)
# plot time at which slab is fully saturated
sat_time_Cu = (Cu_thickness)**2/np.pi/Cu_D
ax.axvline(sat_time_Cu, color='g', linestyle='--', linewidth=1.5, label='Saturation Time')

ax.set_xlabel("Time (s)")
ax.set_ylabel("Length (mum)")
ax.set_xlim(0, 1.2 * sat_time_Cu)
ax.set_ylim(0, np.max(Cu_exact_diff_length))
ax.legend(loc ="best")
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Copper-only: Exact vs Simulated Diffusion Length")
fig.savefig("copper-only-diffusion_length.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Conservation of Mass Check ================
plot_conservation_of_mass(Cu_time, Cu_mass_from_flux, Cu_mass_from_concentration, 0.01, "copper-only-conservation-of-mass.png")

# ==================================== Tungsten only ================================================
# ================ Data Extraction ================
W_time, conc_W, W_exact_diff_length, W_sim_diff_length, W_mass_from_flux, W_mass_from_concentration = read_csv_from_TMAP8("tungsten-only_main.csv", 
                                                                                                                  ["time", "concentration_at_x_W",
                                                                                                                   "exact_diffusion_length", "simulated_diffusion_length",
                                                                                                                   "mass_from_flux", "mass_from_concentration"])
# ================ Plot Concentration over Time ================
fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(W_time, conc_W, color="mediumslateblue")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Concentration (at/μm³)")
ax.set_xlim(0, np.max(W_time))
ax.set_ylim(0, 1.2* np.max(conc_W))
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Tungsten-only concentration vs. time")
fig.savefig("tungsten-only-conc-vs-time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Plot Concentration over Distance ================
W_x, conc_W = read_csv_from_TMAP8(sorted(glob.glob("tungsten-only_vector_W_*.csv"))[-1], ["x", "C_M_W"])

fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(W_x, conc_W, color="darkred")
ax.set_xlabel("Distance (μm)")
ax.set_ylabel("Concentration (at/μm³)")
ax.set_xlim(0, np.max(W_x))
ax.set_ylim(0, 1.2*np.max(conc_W))
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Tungsten-only concentration vs. distance")
fig.savefig("tungsten-only-conc-vs-distance.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Diffusion Front Check ================
# Constant-pressure run: used only for the diffusion-front verification, whose
# analytical solution ell(t) = sqrt(pi*D*t) assumes a step Dirichlet BC.
fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(W_time, W_exact_diff_length, label = "Analytical", color="darkviolet")
ax.plot(W_time, W_sim_diff_length, label = "TMAP8", color="darkred")
annotate_rmspe(
    W_sim_diff_length,
    W_exact_diff_length,
    W_time[-1] / 2,
    W_exact_diff_length[-1] / 4,
)

# plot time at which slab is fully saturated
sat_time_W = (W_thickness)**2/np.pi/W_D
ax.axvline(sat_time_W, color='g', linestyle='--', linewidth=1.5, label='Saturation Time')

ax.set_xlabel("Time (s)")
ax.set_ylabel("Length (mum)")
ax.set_xlim(0, 1.2*sat_time_W) #np.max(W_time))
ax.set_ylim(0, np.max(W_exact_diff_length))
ax.legend(loc ="best")
ax.minorticks_on()
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Tungsten-only: Exact vs Simulated Diffusion Length")
fig.savefig("tungsten-only-diffusion_length.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ================ Conservation of Mass Check ================
plot_conservation_of_mass(W_time, W_mass_from_flux, W_mass_from_concentration, 0.001, "tungsten-only-conservation-of-mass.png")