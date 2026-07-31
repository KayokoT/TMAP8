import os
import glob
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

"""
NOTE: I adapted this case from previous files without much time to test it; there may be flaws.
Assumes run01 is the no trapping case
This file creates:
    - conservation of mass plots for each W_Cu_W_run**
    - penalty verification plots for every run
    - individual plots for H concentration in each material over time
    - individual plots for H concentration in each material over distance
    - a complete distance plot showing concentration across the entire system
    - permeation rate on the downstream graph
"""

# === Input Parameters ======
Cu_thickness = 100  # mum
W_thickness = 25  # mum, thickness of EACH W layer (symmetric stack)
t0 = 10
legend = ["no traps"] 
colors = ["black", "#E69F00", "#0072B2", "#009E73", "red"]
linestyles = ['-', (2, (1, 2.5)), (0, (1, 2.5)), '-', (0, (1, 2.5))]
linewidths = [1.5, 3, 3, 1.5, 3]
# =================================== Functions ===============================
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
    plt.text(x_pos, y_pos, "RMSPE = %.2f %%" % RMSPE, fontweight="bold", fontsize=13)


def plot_conservation_of_mass(t, flux, mass, filename):
    """Plot the absolute percent difference between accumulated boundary flux
    and total mass as a single curve, verifying mass conservation. Data prior
    to t = 1 day is excluded so the metric is not dominated by early-time
    relative noise from the time-integrated flux postprocessor -- analogous to
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
    ax.plot(t, percent_diff, c="r")
    ax.set_xlabel("Time (s)", fontsize=15)
    ax.set_ylabel("Percent difference (%)", fontsize=15)
    ax.set_title("Conservation of Mass Check", fontsize=15)
    ax.set_xlim(left=t[0], right=t[-1])
    ax.set_ylim(bottom=0, top=1)
    plt.grid(which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=13)
    ax.tick_params(axis="both", which="minor", labelsize=12)
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig)


def add_shading(ax, shading):
    """Shade appropriate regions on the concentration vs distance plots

    Args:
        ax (matplotlib.axes.Axes): the axes object to draw on
        shading (list of tuples):
            'label' (str): legend label
            'x0' (float): left x boundary
            'x1' (float): right x boundary
            'color' (str): color
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


def format_ax(ax, xlim):
    """Format concentration vs distance plots"""
    ax.set_xlabel("Distance ($\\mu$m)", fontsize=15)
    ax.set_ylabel(r"Concentration (at/$\mu$m$^3$)", fontsize=15)
    ax.minorticks_on()
    ax.set_xlim(left=xlim[0], right=xlim[1])
    ax.grid(True, which="major", linestyle="--", color="0.65", alpha=0.3)
    ax.tick_params(axis="both", which="major", labelsize=13)
    ax.tick_params(axis="both", which="minor", labelsize=12)
    # ax.legend(loc="best", fontsize=13)


def extend_to_steady_state(time, values, target_end_time, n_extra_points=200):
    """Extend a (time, values) series that has already reached steady state
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


def plot_penalty_check(case, run_id, var_key, gold_key, interface_label, filename_prefix):
    """Plots concentration-ratio vs solubility-ratio over time for a single interface,
    to verify the penalty method is correctly enforcing the solubility jump there.

    Args:
        case (dict): dict with 'time' and the ratio arrays keyed by var_key/gold_key
        run_id (str): identifies this run (e.g. "W_Cu_W_run01_main") -- used in title/filename
        var_key (str): key into `case` for the concentration-ratio time series
        gold_key (str): key into `case` for the solubility-ratio time series
        interface_label (str): e.g. "W1/Cu" -- used in title
        filename_prefix (str): e.g. "W_Cu_W-penalty"
    """
    var_sol = case[var_key]
    gold_sol = case[gold_key]
    time = case["time"]

    fig, ax = plt.subplots(figsize=[6.5, 5.5])
    ax.plot(time, var_sol, label=f"Concentration ({run_id})", linestyle="-")
    ax.plot(time, gold_sol, label=f"Solubility ({run_id})", linestyle="--")

    idx_candidates = np.where(time >= t0)[0]
    if len(idx_candidates) == 0:
        print(f"Warning: no time points >= t0={t0} for {interface_label} ({run_id}); skipping annotation.")
    else:
        idx = idx_candidates[0]
        annotate_rmspe(
            var_sol[idx:],
            gold_sol[idx:],
            0.5 * time[-1],
            1.05 * gold_sol[1]
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Ratio (-)")
    ax.set_ylim(bottom=0, top=1.2 * gold_sol[-1])
    ax.set_xlim(left=t0)
    ax.set_title(f"Comparison of solubility and concentration ratios ({interface_label}, {run_id})")
    ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.legend(loc="best")

    plt.savefig(f"{filename_prefix}-{run_id}-check.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


# ====== Load input files ======

# Main files: concentration/flux/mass time series
main_csvs = sorted(
    glob.glob("W_Cu_W_run[0-9][0-9]_main.csv"),
    key=lambda x: int(re.search(r'run(\d+)', x).group(1))
)

conc_vs_t = {}
for m in main_csvs:
    run = os.path.splitext(os.path.basename(m))[0]  # run01, run02, etc.
    (time, conc_W1, conc_Cu, conc_W2,
     gold_sol_1, var_sol_1, gold_sol_2, var_sol_2,
     right_outflux, mass_from_con, mass_from_flux) = read_csv_from_TMAP8(
        m,
        ["time", "concentration_at_x_W1", "concentration_at_x_Cu", "concentration_at_x_W2",
         "gold_solubility_ratio_1", "variable_ratio_1",
         "gold_solubility_ratio_2", "variable_ratio_2",
         "right_outflux", "mass_from_concentration", "mass_from_flux"]
    )

    conc_vs_t[run] = {
        "time": time,
        "conc_W1": conc_W1,
        "conc_Cu": conc_Cu,
        "conc_W2": conc_W2,
        "sol_ratio_1": gold_sol_1,
        "conc_ratio_1": var_sol_1,
        "sol_ratio_2": gold_sol_2,
        "conc_ratio_2": var_sol_2,
        "mass_from_con": mass_from_con,
        "mass_from_flux": mass_from_flux,
        "right_outflux": right_outflux,
    }

# Files for plotting concentration over distance (most recently modified per run/layer)
dist_csvs = glob.glob("W_Cu_W_run*_*_*.csv")
dist_pattern = re.compile(r"W_Cu_W_(run\d+)_(W1|Cu|W2)_\d+\.csv$")
var_map = {"W1": "C_M_W1", "Cu": "C_M_Cu", "W2": "C_M_W2"}

best_dist_file = {}
for f in dist_csvs:
    m = dist_pattern.search(f)
    if not m:
        continue
    run, layer = m.group(1), m.group(2)
    mtime = os.path.getmtime(f)
    key = (run, layer)
    if key not in best_dist_file or mtime > best_dist_file[key][0]:
        best_dist_file[key] = (mtime, f)

conc_vs_x = {}
for (run, layer), (_, f) in best_dist_file.items():
    var = var_map[layer]
    x, conc = read_csv_from_TMAP8(f, ["x", var])
    conc_vs_x[(run, layer)] = (x, conc)

# sorted list of runs available for distance plots (used by the multi-run plots below)
runs_x = sorted(
    set(r for (r, layer) in conc_vs_x.keys()),
    key=lambda s: int(re.search(r'\d+', s).group())
)

# ========= Concentration as a function of time in first W layer (W1) ===================
all_time = [data["time"] for data in conc_vs_t.values()]
t_max = np.max(np.concatenate(all_time))

all_conc_W1 = [data["conc_W1"] for data in conc_vs_t.values()]
all_conc_Cu = [data["conc_Cu"] for data in conc_vs_t.values()]
all_conc_W2 = [data["conc_W2"] for data in conc_vs_t.values()]
c_max_W1 = np.max(np.concatenate(all_conc_W1))
c_max_Cu = np.max(np.concatenate(all_conc_Cu))
c_max_W2 = np.max(np.concatenate(all_conc_W2))

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

for i, (run, data) in enumerate(conc_vs_t.items()):
    t_ext, conc_ext = extend_to_steady_state(data["time"], data["conc_W1"], t_max)
    ax.plot(t_ext, conc_ext, label=legend[i % len(legend)], linestyle=linestyles[i % len(linestyles)], color=colors[i % len(colors)])

ax.set_xlabel("Time (s)", fontsize=15)
ax.set_ylabel(r"Concentration (at/$\mu$m$^3$)", fontsize=15)
ax.legend(loc="best", fontsize=13)
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)
ax.set_title("Concentration of H in W1 Layer over Time", fontsize=15)
ax.set_xlim(left=0, right=0.0001)
ax.set_ylim(bottom=0, top=1.2 * c_max_W1)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig("W_Cu_W-w1-vs-time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ========= Concentration as a function of time in Cu layer, zoomed to steady state ===================
t_zoom = 0.8 * t_max  # <-- adjust this to wherever steady state actually begins

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

for i, (run, data) in enumerate(conc_vs_t.items()):
    t_ext, conc_ext = extend_to_steady_state(data["time"], data["conc_Cu"], t_max)
    t_ext = np.asarray(t_ext)
    conc_ext = np.asarray(conc_ext)

    mask = t_ext >= t_zoom
    ax.plot(t_ext[mask], conc_ext[mask], label=legend[i % len(legend)], linestyle=linestyles[i % len(linestyles)], color=colors[i % len(colors)])

ax.set_xlabel("Time (s)", fontsize=15)
ax.set_ylabel(r"Concentration (at/$\mu$m$^3$)", fontsize=15)
ax.legend(loc="best", fontsize=13)
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)
ax.set_title("Concentration of H in Cu layer over Time (Steady State)", fontsize=15)
ax.set_xlim(left=t_zoom, right=t_max)
ax.set_ylim(bottom=(1 - 1e-6) * c_max_Cu, top=(1 + 1e-6) * c_max_Cu)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig("W_Cu_W-copper-vs-time-steadystate.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ========= Concentration as a function of time in second W layer (W2) ===================
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

for i, (run, data) in enumerate(conc_vs_t.items()):
    t_ext, conc_ext = extend_to_steady_state(data["time"], data["conc_W2"], t_max)
    ax.plot(t_ext, conc_ext, label=legend[i % len(legend)], linestyle=linestyles[i % len(linestyles)], color=colors[i % len(colors)])

ax.set_xlabel("Time (s)", fontsize=15)
ax.set_ylabel(r"Concentration (at/$\mu$m$^3$)", fontsize=15)
ax.legend(loc="best", fontsize=13)
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)
ax.set_title("Concentration of H in W2 Layer over Time", fontsize=15)
ax.set_xlim(left=0, right=t_max)
ax.set_ylim(bottom=0, top=1.2 * c_max_W2)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
plt.savefig("W_Cu_W-w2-vs-time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Concentration as a function of distance ============
# region boundaries
x_W1Cu = W_thickness
x_CuW2 = W_thickness + Cu_thickness
x_end = 2 * W_thickness + Cu_thickness

regions = {
    "W1": lambda x: (x >= 0) & (x < x_W1Cu),
    "Cu": lambda x: (x >= x_W1Cu) & (x < x_CuW2),
    "W2": lambda x: (x >= x_CuW2) & (x <= x_end),
}

shading = [
    ("Tungsten", 0, x_W1Cu, "gray"),
    ("Copper", x_W1Cu, x_CuW2, "orange"),
    ("Tungsten", x_CuW2, x_end, "gray"),
]

# --- Combined plot: full system, all runs overlaid ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
for i, run in enumerate(runs_x):
    if i == 0:
        x_W1_, c_W1_ = conc_vs_x[(run, "W1")]
        x_Cu_, c_Cu_ = conc_vs_x[(run, "Cu")]
        x_W2_, c_W2_ = conc_vs_x[(run, "W2")]

        # W1, then Cu, then W2
        x_full = np.concatenate([x_W1_, x_Cu_, x_W2_])
        c_full = np.concatenate([c_W1_, c_Cu_, c_W2_])

        ax.plot(
            x_full, c_full,
            linestyle=linestyles[i % len(linestyles)],
            color=colors[i % len(colors)],
        )

for x_interface in (x_W1Cu, x_CuW2):
    ax.axvline(x_interface, color='r', linestyle='--', linewidth=1.5)
    ax.text(
        x_interface - 2, 0.75, 'W–Cu interface',
        transform=ax.get_xaxis_transform(),
        ha='center', va='bottom',
        fontsize=10, color='r',
        rotation=90, rotation_mode='anchor'
    )

add_shading(ax, shading)
format_ax(ax, (0, x_end))
ax.set_yscale("log")
ax.set_title("H Concentration over Distance", fontsize=15)
plt.savefig("W_Cu_W-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- First W layer only (W1) ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
for i, run in enumerate(runs_x):
    x, c = conc_vs_x[(run, "W1")]
    mask = regions["W1"](x)
    ax.plot(
        x[mask], c[mask],
        label=legend[i % len(legend)],
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
    )

add_shading(ax, [shading[0]])
format_ax(ax, (0, x_W1Cu))
ax.legend(loc="best")
ax.set_title("H Concentration over Distance", fontsize=15)
plt.savefig("W_Cu_W-w1-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- Copper layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
for i, run in enumerate(runs_x):
    x, c = conc_vs_x[(run, "Cu")]
    mask = regions["Cu"](x)
    ax.plot(
        x[mask], c[mask],
        label=legend[i % len(legend)],
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
    )

add_shading(ax, [shading[1]])
format_ax(ax, (x_W1Cu, x_CuW2))
ax.legend(loc="best")
ax.set_title("H Concentration over Distance", fontsize=15)
plt.savefig("W_Cu_W-copper-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- Second W layer only (W2) ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
for i, run in enumerate(runs_x):
    x, c = conc_vs_x[(run, "W2")]
    mask = regions["W2"](x)
    ax.plot(
        x[mask], c[mask],
        label=legend[i % len(legend)],
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
    )

add_shading(ax, [shading[2]])
format_ax(ax, (x_CuW2, x_end))
ax.legend(loc="best")
ax.set_title("H Concentration over Distance", fontsize=15)
plt.savefig("W_Cu_W-w2-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Verify Penalty ============
# Plot solubility ratio and concentration ratio to ensure the jump is properly enforced
# at BOTH interfaces (W1/Cu and Cu/W2)
for i, (run_id, case) in enumerate(conc_vs_t.items()):
    label = legend[i % len(legend)]
    plot_penalty_check(
        case,
        run_id=label,
        var_key="conc_ratio_1",
        gold_key="sol_ratio_1",
        interface_label="W1/Cu",
        filename_prefix="W_Cu_W-penalty-interface1",
    )
    plot_penalty_check(
        case,
        run_id=label,
        var_key="conc_ratio_2",
        gold_key="sol_ratio_2",
        interface_label="Cu/W2",
        filename_prefix="W_Cu_W-penalty-interface2",
    )

# ============ Verify Conservation of Mass ============
for run, data in conc_vs_t.items():
    title = f"W_Cu_W-cons-mass_{run}"
    plot_conservation_of_mass(
        data["time"],
        data["mass_from_flux"],
        data["mass_from_con"],
        f"{title}.png"
    )

# ============ Permeation flux ============
# Compares the permeation rate between trapping and no trapping
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

# Reference case (no trapping)
ref_run = list(conc_vs_t.keys())[0]
t_ref = conc_vs_t[ref_run]["time"]
perm_ref = conc_vs_t[ref_run]["right_outflux"]

for i, (run, data) in enumerate(conc_vs_t.items()):
    t_ext, perm_ext = extend_to_steady_state(data["time"], data["right_outflux"], t_max)

    if run == ref_run:
        label = f"{legend[i]} (RMSPE = 0.00%)"
    else:
        # Interpolate onto the reference time grid
        perm_interp = np.interp(t_ref, t_ext, perm_ext)
        rmspe = compute_rmspe(perm_interp, perm_ref)
        label = f"{legend[i]} (RMSPE = {rmspe:.2f}%)"

    ax.plot(
        t_ext, perm_ext,
        label=label,
        linestyle=linestyles[i % len(linestyles)],
        color=colors[i % len(colors)],
        linewidth=linewidths[i % len(linewidths)]
    )

ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
ax.set_xlabel("Time", fontsize=15)
ax.set_ylabel("Permeation Rate (H atoms/$\\mu$m$^2$/s)", fontsize=15)
ax.legend(loc="lower right", fontsize=10.5)
ax.tick_params(axis="both", which="major", labelsize=11)
ax.tick_params(axis="both", which="minor", labelsize=11)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.set_title("Permeation over Time", fontsize=15)
ax.minorticks_on()
plt.savefig("W_Cu_W-permeation-rate.png", bbox_inches="tight", dpi=300)
plt.close(fig)
