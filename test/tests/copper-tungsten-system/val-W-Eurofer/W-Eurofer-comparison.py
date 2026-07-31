from ast import main
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
 This file compares the TMAP8 simulation results for the W-Eurofer system (W-Eurofer.i) at three pressures
 (107.6kPa, 101.4kPa, 95.3 kPa, changed in the W-Eurofer.params file) to the experimental results from Zajec et al's paper, 
 "Hydrogen diffusive transport parameters in W coating for fusion applications", Journal of Nuclear Materials,
 http://dx.doi.org/10.1016/j.jnucmat.2011.02.039.
 
 Specifically, this file compares the permeation rate of the simulations to the experiments.
 It was concluded that there are large discrepancies, likely due to the ambiguous derivation of the
 solubilities and diffusivities in the paper, as well as the fact that surface kinetics and traps
 confound the pressure-dependent results. Thus, this validation effort was abandoned.
 
 In addition to validation, this file outputs conservation of mass and penalty checks on the simulation.
 """
# =================================== Function ===============================
def read_csv_from_TMAP8(parameter_names, csv_folder):
    """Read simulation output columns from a CSV file into a numpy array
    From AUNF Mini-Canister Example Case

    Args:
        parameter_names (list of str): column names to extract, in desired order
        csv_folder (str) path to data

    Returns:
        ndarray: 2D array of shape (len(parameter_names), n_timesteps)
    """
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


def plot_conservation_of_mass(t, flux, mass, filename):
    """Plot the absolute percent difference between accumulated boundary flux
    and total mass as a single curve, verifying mass conservation. Data prior
    to t = 1 day is excluded so the metric is not dominated by early-time
    relative noise from the time-integrated flux postprocessor — analogous to
    the early-time slicing convention used for RMSPE in the ver-* cases.

    Args:
        t (float, ndarray): time array in days
        flux (float, ndarray): accumulated boundary flux in µmol H
        mass (float, ndarray): total H mass in domain in µmol H
        filename (str): output PNG filename
    """
    t = np.asarray(t)
    flux = np.asarray(flux, dtype=float)
    mass = np.asarray(mass, dtype=float)
    idx = np.where(t > 0.001)[0][0]
    t = t[idx:]
    percent_diff = np.abs(flux[idx:] / mass[idx:] - 1.0) * 100.0

    fig = plt.figure(figsize=[6.5, 5.5])
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(t, percent_diff, c="tab:gray")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Percent difference (%)")
    ax.set_xlim(left=t[0], right=t[-1])
    ax.set_ylim(bottom=0)
    plt.title("Conservation of Mass Check")
    plt.grid(which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig)

def add_shading(ax, shading):
    """Shade appropriate regions on the concentration vs distance plots
    Args:
        ax (matplotlib.axes.Axes): the axes object to draw on
        shading (list of tuples): 
            'label'    (str): legend label
            'x0'    (float): left x boundary
            'x1'      (float):   right x boundary
            'color'   (str):     color
        filename (str): output PNG filename
    """
    for label, x0, x1, color in shading:
        ax.axvspan(x0, x1, color=color, alpha=0.15)
        mid = 0.5 * (x0 + x1)
        ax.text(
            mid, 0.5, label,
            transform=ax.get_xaxis_transform(),
            ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="dimgray",
        )


def extend_to_steady_state(time, values, target_end_time, n_extra_points=200):
    """
    Extend a (time, values) series that has already reached steady state
    by holding the final value constant out to target_end_time.
    """
    time = np.asarray(time)
    values = np.asarray(values)

    if time[-1] >= target_end_time:
        return time, values

    extra_time = np.linspace(time[-1], target_end_time, n_extra_points)[1:]  # skip duplicate endpoint
    extra_values = np.full_like(extra_time, values[-1])

    time_extended = np.concatenate([time, extra_time])
    values_extended = np.concatenate([values, extra_values])

    return time_extended, values_extended

def numerical_solution_on_experiment_input(experiment_input, tmap_input, tmap_output):
    """Get new numerical solution based on the experimental input data points

    Args:
        experiment_input (float, ndarray): experimental input data points
        tmap_input (float, ndarray): numerical input data points
        tmap_output (float, ndarray): numerical output data points

    Returns:
        float, ndarray: updated tmap_output based on the data points in experiment_input
    """
    new_tmap_output = np.zeros(len(experiment_input))
    for i in range(len(experiment_input)):
        left_limit = np.argwhere((np.diff(tmap_input < experiment_input[i])))[0][0]
        right_limit = left_limit + 1
        new_tmap_output[i] = (experiment_input[i] - tmap_input[left_limit]) / (
            tmap_input[right_limit] - tmap_input[left_limit]
        ) * (tmap_output[right_limit] - tmap_output[left_limit]) + tmap_output[
            left_limit
        ]
    return new_tmap_output

def format_ax(ax, xlim):
    """Format concentration vs distance plots
    Args:
        ax (matplotlib.axes.Axes): the axes object to draw on
        xlim (tuple):
            left (float): left x-limit
            right (float): right x-limit
    """
    ax.set_xlabel("Distance (mum)")
    ax.set_ylabel(r"Concentration (at/mum$^3$)")
    ax.minorticks_on()
    ax.set_xlim(left=xlim[0], right=xlim[1])
    ax.grid(True, which="major", linestyle="--", color="0.65", alpha=0.3)


# ====== Load experimental + simulation data for each pressure ======
pressures = ["107600", "101400", "95300"]
sim_parameter_names = [
    "time","gold_solubility_ratio", "variable_ratio", "right_outflux_H2",
    "mass_from_concentration", "mass_from_flux", "concentration_at_x_W", "concentration_at_x_Euro"
]
cases = []
for p in pressures:
    exp_time, exp_perm = read_csv_from_TMAP8(
        [f"{p}Pa-time (s)", f"{p}Pa-j (molH2/m^2/s)"], "./gold/experiment_data_permeation_transient.csv"
    )

    mask = ~np.isnan(exp_perm)
    exp_time = exp_time[mask]
    exp_perm = exp_perm[mask]

    (time, gold_sol, var_sol,
     permeation, mass_from_con, mass_from_flux,
     concentration_at_x_W, concentration_at_x_Euro) = read_csv_from_TMAP8(sim_parameter_names, f"./W_Eurofer_{p}_main.csv"
    )

    cases.append(dict(
        time=time,
        var_sol=var_sol,
        gold_sol=gold_sol,
        pressure=p,
        permeation=permeation,
        experiment_time=exp_time,
        experiment_perm=exp_perm,
        mass_from_flux = mass_from_flux,
        mass_from_con = mass_from_con,
        concentration_at_x_W = concentration_at_x_W,
        concentration_at_x_Euro=concentration_at_x_Euro
    ))
# ==================================== Permeation flux ===============================================
fig, ax = plt.subplots(figsize=[6.5, 5.5])

for i, case in enumerate(cases):
    line, = ax.plot(
        case["time"], case["permeation"],
        label=f"TMAP8 {case['pressure']} Pa", linestyle="-",
    )
    ax.plot(
        case["experiment_time"], case["experiment_perm"],
        label=f"Experiment {case['pressure']} Pa", linestyle="--",
        color=line.get_color(),
    )

    tmap_flux_for_rmspe = numerical_solution_on_experiment_input(
        case["experiment_time"], case["time"], case["permeation"]
    )
    RMSPE = compute_rmspe(tmap_flux_for_rmspe, case["experiment_perm"])

    ax.text(
        1e4 / 3600.0, (40 - 5 * i) * 1e15,
        "RMSPE = %.2f %%" % RMSPE,
        fontweight="bold", color=line.get_color(),
    )

ax.set_xlabel("Time (s)")
ax.set_ylabel("Permeation rate (mol H$_2$ m$^{-2}$ s$^{-1}$)")
ax.set_title("Permeation rate comparison")
ax.legend(loc="best")
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()

fig.savefig("W-Eurofer-permeation-all-pressures.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Verify Penalty ============
# Plot solubility ratio and concentration ratio to ensure the jump is properly enforced
t0 = 1000 # choose where concentration stabilizes
# Plot solubility ratio and concentration ratio to ensure the jump is properly enforced

for case in cases:
    fig = plt.figure(figsize=[6.5, 5.5])
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(case["time"], case["var_sol"], label="Concentration ratio", color="tab:red", linestyle="-")
    ax.plot(case["time"], case["gold_sol"],  label="Solubility ratio", color="k", linestyle="--")
    idx = np.where(time >= t0)[0][0]
    annotate_rmspe(case["var_sol"][idx:], case["gold_sol"][idx:], 0.01, 1.05*case["gold_sol"][1])
    #ax.set_ylim(bottom=0, top=1.2 * case["gold_sol"][1])
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Ratio (-)")
    ax.legend(loc="best")
    plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.set_title("Comparison of solubility and concentration ratios (W/Euro)")
    ax.minorticks_on()
    plt.savefig(f"W-Eurofer-penalty-check-{case["pressure"]}.png", bbox_inches="tight", dpi=300)

# ==================================== Verify Concentration of Mass ================================================
for case in cases:
    plot_conservation_of_mass(
        case["time"], case["mass_from_flux"], case["mass_from_con"],
        f"W-Eurofer-{case['pressure']}-cons-of-mass.png"
    )


# # ============ Comparison of concentration as a function of time ============
fig, ax = plt.subplots(figsize=[6.5, 5.5])

for i, case in enumerate(cases):
    line, = ax.plot(
        case["time"], case["concentration_at_x_W"],
        label=f"{case['pressure']} Pa", linestyle="-",
    )

ax.set_xlabel("Time (s)")
ax.set_ylabel("H Concentration (mol H m$^{-3}$)")
ax.set_title("Concentration in Tungsten over Time")
ax.legend(loc="best")
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()

fig.savefig("W-Eurofer-W-conc-vs-t.png", bbox_inches="tight", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=[6.5, 5.5])

for i, case in enumerate(cases):
    line, = ax.plot(
        case["time"], case["concentration_at_x_Euro"],
        label=f"{case['pressure']} Pa", linestyle="-",
    )

ax.set_xlabel("Time (s)")
ax.set_ylabel("H Concentration (mol H m$^{-3}$)")
ax.set_title("Concentration in Eurofer over Time")
ax.legend(loc="best")
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()

fig.savefig("W-Eurofer-E-conc-vs-t.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# # ============ Comparison of concentration as a function of distance ============
def latest_file(pattern):
    candidates = glob.glob(pattern)
    if not candidates:
        raise FileNotFoundError(f"No files found matching {pattern}")
    return max(candidates, key=os.path.getmtime)

Eurofer_file = latest_file("W_Eurofer_95300_Euro_*.csv")
W_file = latest_file("W_Eurofer_95300_W_*.csv")

x_Eurofer, conc_Eurofer = read_csv_from_TMAP8(["x", "C_M_Euro"], Eurofer_file)
x_W, conc_W = read_csv_from_TMAP8(["x", "C_M_W"], W_file)

# region boundaries
W_thickness = 5 * 1e-6
Euro_thickness = 0.5 * 1e-3
x_WEuro = W_thickness
x_Euro = W_thickness + Euro_thickness

# masking functions
regions = {
    "W":       lambda x: x < x_WEuro,
    "Eurofer": lambda x: (x >= x_WEuro) & (x < x_Euro),
}

shading = [
    ("Tungsten", 0,       x_WEuro, "forestgreen"),
    ("Eurofer",  x_WEuro, x_Euro,  "palegreen"),
]

# --- Combined plot ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
mask_W = regions["W"](x_W)
mask_Eurofer = regions["Eurofer"](x_Eurofer)
ax.plot(x_W[mask_W], conc_W[mask_W], c="k")
ax.plot(x_Eurofer[mask_Eurofer], conc_Eurofer[mask_Eurofer], c="k")
ax.axvline(x_WEuro, color='g', linestyle='--', linewidth=1.5, label='W-Eurofer interface')
add_shading(ax, shading)
format_ax(ax, (0, x_Euro))
ax.set_title("Concentration over Distance — Full System")
fig.savefig("W-Eurofer-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- Eurofer layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
mask = regions["Eurofer"](x_Eurofer)
ax.plot(x_Eurofer[mask], conc_Eurofer[mask], c="k")
add_shading(ax, [shading[1]])  # Eurofer shading
format_ax(ax, (x_WEuro, x_Euro))
ax.set_title("Concentration over Distance — Eurofer Layer")
fig.savefig("W-Eurofer-Eurofer-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- W layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
mask = regions["W"](x_W)
ax.plot(x_W[mask], conc_W[mask], c="k")
add_shading(ax, [shading[0]])  # Tungsten shading
format_ax(ax, (0, x_WEuro))
ax.set_title("Concentration over Distance — W Layer")
fig.savefig("W-Eurofer-tungsten-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

