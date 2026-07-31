import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
import os
import glob
# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

"""
 This file compares the TMAP8 simulation results for the CVD W-Cu system (val-W-Cu.i), as well as just CVD-W and copper (val-Cu.i) 
 at various temperatures to the experimental results from Sun et al's paper, 
 "Deuterium plasma-driven permeation in chemical vapor deposition tungsten-copper composite", Nuclear Fusion,
 https://doi.org/10.1088/1741-4326/add171.
 
 Specifically, this file compares the normalized permeation rate of the simulations to the experiments.

 In addition to validation, this file outputs conservation of mass and penalty checks on the simulation.

 To properly run, this file requires all 5 Cu-only temperature cases, all 4 W-only cases, and all 3 W-Cu cases
 to be run.

 bilayer-perm-validation-741.png refers to the comparison of the simulated data against the actual 
 plotted permeation (D/m^2/s), whereas bilayer-permeation-741.png refers to the comparison of the simulated data 
 against the normalized 741 case multiplied by its steady state; thus small discrepancies between the two appear 
 due to data digitization uncertainties
 """
# =================================== FUNCTIONS ===============================
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

def compute_steady_state_rmspe(simulated, reference):
    """Compute the Root Mean Square Percentage Error between two arrays,
    using only the last quarter of each array (i.e., once the simulation
    or reference series is more than 75% complete). This focuses the
    metric on steady-state agreement rather than transient behavior.

    Args:
        simulated (ndarray): simulated values
        reference (ndarray): reference values used as the denominator

    Returns:
        float: RMSPE in percent, computed over the last quarter of the data
    """
    simulated = np.asarray(simulated)
    reference = np.asarray(reference)

    n = len(reference)
    start = (3 * n) // 4

    sim_quarter = simulated[start:]
    ref_quarter = reference[start:]

    RMSE = np.sqrt(np.mean((sim_quarter - ref_quarter) ** 2))
    return RMSE * 100 / np.mean(ref_quarter)

def annotate_rmspe(simulated, reference, x_pos, y_pos):
    """Compute RMSPE and annotate it as bold text on the current matplotlib axes

    Args:
        simulated (float, ndarray): simulated values
        reference (float, ndarray): reference values used as the denominator
        x_pos (float): x-coordinate of the annotation
        y_pos (float): y-coordinate of the annotation
    """
    RMSPE = compute_rmspe(simulated, reference)
    plt.text(x_pos, y_pos, "RMSPE = %.2f %%" % RMSPE, fontweight="bold", color='r', fontsize = 13)


def plot_conservation_of_mass(t, flux, mass, filename):
    """Plot the absolute percent difference between accumulated boundary flux
    and total mass as a single curve, verifying mass conservation. Data prior
    to t = 100 s is excluded so the metric is not dominated by early-time
    relative noise from the time-integrated flux postprocessor — analogous to
    the early-time slicing convention used for RMSPE in the ver-* cases.

    Args:
        t (float, ndarray): time array in seconds
        flux (float, ndarray): accumulated boundary flux in µmol H
        mass (float, ndarray): total H mass in domain in µmol H
        filename (str): output PNG filename
    """
    t = np.asarray(t)
    flux = np.asarray(flux, dtype=float)
    mass = np.asarray(mass, dtype=float)
    idx = np.where(t > 100)[0][0]
    t = t[idx:]
    percent_diff = np.abs(flux[idx:] / mass[idx:] - 1.0) * 100.0

    fig = plt.figure(figsize=[6.5, 5.5])
    gs = gridspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0])
    ax.plot(t, percent_diff, c="tab:gray")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Percent difference (%)")
    ax.set_xlim(left=t[0], right=t[-1])
    ax.set_ylim(bottom=0, top = 2)
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


def plot_penalty_check(case, var_key, gold_key, interface_label, filename_prefix):
    """
    Plots concentration-ratio vs solubility-ratio over time for a single interface,
    to verify the penalty method is correctly enforcing the solubility jump there.

    case: dict with 'time', 'temperature', and the ratio arrays keyed by var_key/gold_key
    var_key, gold_key: dict keys into `case` for the concentration-ratio and
                        solubility-ratio time series for this specific interface
    interface_label: e.g. "W/CuEntry" or "CuEntry/Cu" -- used in title
    filename_prefix: e.g. "bilayer-penalty-W-CuEntry"
    """
    var_sol = case[var_key]
    gold_sol = case[gold_key]
    time = case["time"]
    T = case["temperature"]

    fig, ax = plt.subplots(figsize=[6.5, 5.5])
    ax.plot(time, var_sol, label=f"Concentration ({T} K)", linestyle="-")
    ax.plot(time, gold_sol, label=f"Solubility ({T} K)", linestyle="--")

    idx_candidates = np.where(time >= t0)[0]
    if len(idx_candidates) == 0:
        print(f"Warning: no time points >= t0={t0} for {interface_label} at {T} K; skipping annotation.")
    else:
        idx = idx_candidates[0]
        annotate_rmspe(
            var_sol[idx:],
            gold_sol[idx:],
            0.01,
            1.05 * gold_sol[1]
        )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Ratio (-)")
    ax.set_ylim(bottom=0, top=1.2 * gold_sol[-1])
    ax.set_xlim(left=10)
    ax.set_title(f"Comparison of solubility and concentration ratios ({interface_label})")
    ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.legend(loc="best")

    plt.savefig(f"{filename_prefix}-{T}-check.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
# ===================================== CU SYSTEM ===============================================
# Load experimental + simulation data for each temperature
temperatures_Cu = ["587", "620", "700", "742"] #K
MAIN_COLUMNS = [
    "time", "right_outflux",
    "mass_from_concentration", "mass_from_flux"
]
cases_Cu = []
for t in temperatures_Cu:
    exp_time, exp_perm = read_csv_from_TMAP8(
        [f"time_{t}", f"perm_{t}"], "./gold/Cu-experiment-normperm.csv"
    )

    (time, permeation, mass_from_con, mass_from_flux
    ) = read_csv_from_TMAP8(MAIN_COLUMNS, f"./val_Cu_{t}_main.csv"
    )

    cases_Cu.append(dict(
        time=time,
        temperature=t,
        permeation=permeation,
        experiment_time=exp_time,
        experiment_perm=exp_perm,
        mass_from_flux = mass_from_flux,
        mass_from_con = mass_from_con
    ))

# compare simulated and experimental Cu normalized permeation flux
fig, ax = plt.subplots(figsize=[6.5, 5.5])
for i, case in enumerate(cases_Cu):
    normalized_perm = case["permeation"]/case["permeation"][-1]
    line, = ax.plot(
        case["time"], normalized_perm,
        label=f"TMAP8 {case['temperature']} K", linestyle="-"
    )
    ax.plot(
        case["experiment_time"], case["experiment_perm"],
        label=f"Experiment {case['temperature']} K", linestyle="--",
        color=line.get_color(),
    )

ax.set_xlabel("Time (s)")
ax.set_ylabel("Normalized Permeation (-)")
ax.set_title("Copper Only - Permeation Rate Comparison")
ax.legend(loc="best")
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
fig.savefig("Cu-permeation.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# verify concentration of mass
for case in cases_Cu:
    plot_conservation_of_mass(
        case["time"], case["mass_from_flux"], case["mass_from_con"],
        f"Cu-{case['temperature']}-cons-of-mass.png"
    )
# directly compare Cu simulation against permeation rate vs time (Fig 7) at 741 K
(time_Cu_741, permeation_Cu_741, mass_from_con_Cu_741, mass_from_flux_Cu_741) = read_csv_from_TMAP8(MAIN_COLUMNS, "./val_Cu_741_main.csv"
    )
permeation_Cu_741 = permeation_Cu_741 * 1e12 #convert at/mum^2/s -> at/m^2/s
exp_time_Cu_741, exp_perm_Cu_741 = read_csv_from_TMAP8(
    ["time", "flux"], "./gold/Cu-experiment-flux.csv"
)

fig, ax = plt.subplots(figsize=[6.5, 5.5])

ax.plot(
    time_Cu_741,
    permeation_Cu_741,
    label='simulation',
    linestyle="-"
)
ax.plot(
    exp_time_Cu_741, 
    exp_perm_Cu_741,
    label="experiment", 
    linestyle="--"
    )

tmap_flux_for_rmspe = np.interp(
       exp_time_Cu_741, time_Cu_741, permeation_Cu_741
)

annotate_rmspe(tmap_flux_for_rmspe, exp_perm_Cu_741, 0.5*max(exp_time_Cu_741), 1.2*exp_perm_Cu_741[-1])
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Permeation flux (D/m$^2$/s)")
ax.set_yscale("log")
ax.set_ylim(bottom = 1e14, top = 1e20)
ax.set_xlim(left = 0, right = 250)
ax.set_title("Cu - TMAP8 vs Experimental Permeation Curves", fontsize=15)
ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
ax.legend(loc="best")
plt.savefig("Cu-perm-validation.png", bbox_inches="tight", dpi=300)
plt.close()

# check conservation of mass
plot_conservation_of_mass(
        time_Cu_741, mass_from_flux_Cu_741, mass_from_con_Cu_741,
        "Cu-741-cons-of-mass.png"
    )

# ===================================== W SYSTEM ===============================================
# Load experimental + simulation data for each temperature
temperatures_W = ["643", "742", "790"] #K
cases_W = []
for t in temperatures_W:
    exp_time, exp_perm = read_csv_from_TMAP8(
        [f"time_{t}", f"perm_{t}"], "./gold/CVD-W-experiment-normperm.csv"
    )

    exp_time = np.asarray(exp_time, dtype=float)
    exp_perm = np.asarray(exp_perm, dtype=float)

    mask = ~np.isnan(exp_perm)
    exp_time = exp_time[mask]
    exp_perm = exp_perm[mask]

    (
        time,
        permeation,
        mass_from_con,
        mass_from_flux
    ) = read_csv_from_TMAP8(
        MAIN_COLUMNS,
        f"./val_W_{t}_main.csv"
    )

    cases_W.append(dict(
        time=time,
        temperature=t,
        permeation=permeation,
        experiment_time=exp_time,
        experiment_perm=exp_perm,
        mass_from_flux = mass_from_flux,
        mass_from_con = mass_from_con
    ))


# compare simulated and experimental W normalized permeation flux 
fig, ax = plt.subplots(figsize=[6.5, 5.5])
for i, case in enumerate(cases_W):
    normalized_perm = case["permeation"] / case["permeation"][-1]

    tmap_perm_for_rmspe = np.interp(
        case["experiment_time"],
        case["time"],
        normalized_perm,
    )

    mask = case["experiment_perm"] > 0.05

    RMSPE = compute_rmspe(
    tmap_perm_for_rmspe[mask],
    case["experiment_perm"][mask],
)
    line, = ax.plot(
        case["time"],
        normalized_perm,
        label=f"TMAP8 {case['temperature']} K (RMSPE={100*RMSPE:.2f}%)",
        linestyle="-",
    )

    ax.plot(
        case["experiment_time"],
        case["experiment_perm"],
        label=f"Experiment {case['temperature']} K",
        linestyle="--",
        color=line.get_color(),
    )

ax.set_xlabel("Time (s)")
ax.set_ylabel("Normalized Permeation")
ax.set_title("Tungsten Only - Permeation Rate Comparison")
ax.legend(loc="best")
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
fig.savefig("W-permeation.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# verify conservation of mass
for case in cases_W:
    plot_conservation_of_mass(
        case["time"], case["mass_from_flux"], case["mass_from_con"],
        f"W-{case['temperature']}-cons-of-mass.png"
    )

# directly compare W simulation against permeation rate vs time (Fig 7) at 741 K
(time_W_741, permeation_W_741, mass_from_con_W_741, mass_from_flux_W_741) = read_csv_from_TMAP8(MAIN_COLUMNS, "./val_W_741_main.csv"
    )
permeation_W_741 = permeation_W_741 * 1e12 #convert at/mum^2/s -> at/m^2/s
#Load experimental data
exp_time_W_741, exp_perm_W_741 = read_csv_from_TMAP8(
        ["time", "flux"], "./gold/CVD-W-experiment-flux.csv"
    )
fig, ax = plt.subplots(figsize=[6.5, 5.5])

ax.plot(
    time_W_741,
    permeation_W_741,
    label="simulation",
    linestyle="-"
)
ax.plot(
    exp_time_W_741, exp_perm_W_741,
    label="experiment", 
    linestyle="--"
    )
tmap_flux_for_rmspe = np.interp(
       exp_time_W_741, time_W_741, permeation_W_741
)

annotate_rmspe(tmap_flux_for_rmspe, exp_perm_W_741,
        1000, 4e17)
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Permeation flux (D/m$^2$/s)")
ax.set_yscale("log")
ax.set_ylim(bottom = 1e14, top = 1e20)
ax.set_xlim(left = 0, right = 4000)
ax.set_title("W - TMAP8 vs Experimental Permeation Curves", fontsize=15)
ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
ax.legend(loc="best")
plt.savefig("W-perm-validation.png", bbox_inches="tight", dpi=300)
plt.close()

# check conservation of mass
plot_conservation_of_mass(
        time_W_741, mass_from_flux_W_741, mass_from_con_W_741,
        "W-741-cons-of-mass.png"
    )
# ===================================== W-CUENTRY-CU SYSTEM ==============================================
# Accounts for an inner region with different solubility between pure W and pure Cu where Cu enters W
# Load experimental + simulation data for each temperature
temperatures_bilayer = ["660", "700", "741"] #K

BILAYER_MAIN_COLUMNS = [
    "time", "right_outflux",
    "gold_solubility_ratio_1", "variable_ratio_1",   # W / CuEntry interface
    "gold_solubility_ratio_2", "variable_ratio_2",   # CuEntry / Cu interface
    "mass_from_concentration", "mass_from_flux"
]

cases_bilayer = []
for t in temperatures_bilayer:
    exp_time, exp_perm = read_csv_from_TMAP8(
        [f"time_{t}", f"perm_{t}"], "./gold/CVD-W-Cu-experiment-normperm.csv"
    )

    mask = ~np.isnan(exp_perm)
    exp_time = exp_time[mask]
    exp_perm = exp_perm[mask]

    (time, permeation,
     gold_sol_1, var_sol_1,
     gold_sol_2, var_sol_2,
     mass_from_con, mass_from_flux) = read_csv_from_TMAP8(
        BILAYER_MAIN_COLUMNS, f"./val_W_entryCu_Cu_{t}_main.csv"
    )

    cases_bilayer.append(dict(
        time=time,
        temperature=t,
        permeation=permeation,
        experiment_time=exp_time,
        experiment_perm=exp_perm,
        mass_from_flux=mass_from_flux,
        mass_from_con=mass_from_con,
        gold_sol_1=gold_sol_1,
        var_sol_1=var_sol_1,
        gold_sol_2=gold_sol_2,
        var_sol_2=var_sol_2,
    ))
# compare simulated and experimental CVD W/Cu normalized permeation flux
fig, ax = plt.subplots(figsize=[6.5, 5.5])
for i, case in enumerate(cases_bilayer):
    normalized_perm = case["permeation"]/case["permeation"][-1]
    line, = ax.plot(
        case["time"], normalized_perm,
        label=f"TMAP8 {case['temperature']}K", linestyle="-",
    )
    ax.plot(
        case["experiment_time"], case["experiment_perm"],
        label=f"Experiment {case['temperature']}K", linestyle="--",
        color=line.get_color(),
    )

ax.set_xlabel("Time (s)", fontsize = 15)
ax.set_ylabel("Normalized Permeation (-)", fontsize = 15)
ax.set_title("W/Cu - Normalized Permeation Rate Comparison", fontsize = 15)
ax.legend(loc="best", fontsize = 13)
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)
ax.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
fig.savefig("bilayer-normalized-permeation-all-temps.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# compare simulated and experimental CVD W/Cu permeation flux, NOT normalized
# uses Fig 3 for the normalized values, multiplied by the steady state values in Fig 5 to unscale it
perm_data = pd.read_csv("./gold/permflux-vs-temp-experiment.csv")

for case in cases_bilayer:
    temperature = case["temperature"]

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    # TMAP permeation profile
    tmap_time = case["time"]
    tmap_perm = case["permeation"] * 1e12  # D/mum^2/s -> D/m^2/s

    # get normalized experimental permeation profile
    exp_time = case["experiment_time"]
    exp_norm_perm = case["experiment_perm"]

    # Convert normalized profile to absolute flux using this case's steady-state value
    exp_ss = perm_data.loc[
        perm_data["CVD-W-Cu-temp"] == int(temperature),
        "CVD-W-Cu-ss"
    ].values[0]
    exp_perm = exp_norm_perm * exp_ss

    # calculate RMSPE
    tmap_perm_interp = np.interp(
        exp_time,
        tmap_time,
        tmap_perm,
    )
    # plot comparison
    ax.plot(
        tmap_time,
        tmap_perm,
        label="TMAP8",
        linewidth=2,
    )

    ax.plot(
        exp_time,
        exp_perm,
        linestyle="--",
        label="Experiment",
    )
    annotate_rmspe(tmap_perm_interp, exp_perm, 0.5*max(tmap_time), 1.2*tmap_perm[-1])
    print(f"the steady state RMSPE for the bilayer case at {temperature} K is s {compute_steady_state_rmspe(tmap_perm_interp, exp_perm)} %")

    ax.set_xlabel("Time (s)", fontsize=15)
    ax.set_ylabel("Permeation flux (D/m$^2$/s)", fontsize=15)
    ax.set_xlim(left=0, right=max(tmap_time))
    ax.set_ylim(bottom=1e14, top=1e20)
    ax.set_yscale("log")
    ax.set_title(f"Deuterium permeation at {temperature} K", fontsize=15)
    ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.legend(loc="best", fontsize=13)

    plt.savefig(
        f"bilayer-permeation-{temperature}.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close(fig)


colors = ["slategray", "crimson", "royalblue", "darkviolet"]

fig, ax = plt.subplots(figsize=(7.5, 6))

for i, case in enumerate(cases_bilayer):
    temperature = case["temperature"]
    color = colors[i % len(colors)]

    # TMAP permeation profile
    tmap_time = case["time"]
    tmap_perm = case["permeation"] * 1e12  # D/mum^2/s -> D/m^2/s

    # get normalized experimental permeation profile
    exp_time = case["experiment_time"]
    exp_norm_perm = case["experiment_perm"]

    # Convert normalized profile to absolute flux using this case's steady-state value
    exp_ss = perm_data.loc[
        perm_data["CVD-W-Cu-temp"] == int(temperature),
        "CVD-W-Cu-ss"
    ].values[0]
    exp_perm = exp_norm_perm * exp_ss

    # interpolate TMAP onto experimental time points and compute RMSPE for the legend
    tmap_perm_interp = np.interp(
        exp_time,
        tmap_time,
        tmap_perm,
    )
    RMSPE = compute_rmspe(
        tmap_perm_interp,
        exp_perm,
    )

    # simulation curve (solid)
    ax.plot(
        tmap_time,
        tmap_perm,
        label=f"TMAP8 {temperature}K", # (RMSPE={100*RMSPE:.1f}%)",
        linestyle="-",
        color=color,
        linewidth=2,
    )

    # experiment curve (dashed, same color as its simulation counterpart)
    ax.plot(
        exp_time,
        exp_perm,
        label=f"Experiment {temperature}K",
        linestyle="--",
        color=color,
        linewidth=1.5,
    )

ax.set_xlabel("Time (s)", fontsize=15)
ax.set_ylabel("Permeation flux (D/m$^2$/s)", fontsize=15)
ax.set_xlim(left=0, right = 16000)
ax.set_ylim(bottom=1e14, top=1e20)
ax.set_yscale("log")
ax.set_title("Deuterium permeation in CVD-W/Cu", fontsize=15)
ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
ax.legend(loc="best", fontsize=11)
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)

plt.savefig(
    "bilayer-permeation-all-temps-combined.png",
    bbox_inches="tight",
    dpi=300,
)
plt.close(fig)

# verify conservation of mass
for case in cases_bilayer:
    plot_conservation_of_mass(
        case["time"], case["mass_from_flux"], case["mass_from_con"],
        f"bilayer-{case['temperature']}-cons-of-mass.png"
    )



# directly compare W/Cu simulation against permeation rate vs time (Fig 7) at 741 K
(time_741, permeation_741, mass_from_con_741, mass_from_flux_741) = read_csv_from_TMAP8(MAIN_COLUMNS, "./val_W_entryCu_Cu_741_main.csv"
    )

permeation_741 = permeation_741 * 1e12 #convert at/mum^2/s -> at/m^2/s
exp_time_741, exp_perm_741 = read_csv_from_TMAP8(
["time", "flux"], "./gold/CVD-W-Cu-experiment-flux.csv"
)

fig, ax = plt.subplots(figsize=[6.5, 5.5])
ax.plot(
    time_741,
    permeation_741,
    label='simulation',
    linestyle="-"
)
ax.plot(
    exp_time_741, exp_perm_741,
    label="experiment", 
    linestyle="--"
    )

tmap_flux_for_rmspe = np.interp(
       exp_time_741, time_741, permeation_741
)
RMSPE = compute_rmspe(tmap_flux_for_rmspe, exp_perm_741)

ax.text(
        4000, 4e18,
        "RMSPE = %.2f %%" % RMSPE,
        fontweight="bold", color="r", fontsize = 13
)

ax.set_xlabel("Time (s)", fontsize=15)
ax.set_ylabel(r"Permeation flux (D/m$^2$/s)", fontsize=15)
ax.set_yscale("log")
ax.set_ylim(bottom=1e14, top=1e20)
ax.set_xlim(left=0, right=7000)
ax.set_title("Deuterium permeation at 741 K", fontsize=15)
ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
ax.minorticks_on()
ax.tick_params(axis="both", which="major", labelsize=13)
ax.tick_params(axis="both", which="minor", labelsize=12)
ax.legend(loc="best", fontsize=13)
plt.savefig("bilayer-perm-validation-741.png", bbox_inches="tight", dpi=300)
plt.close()

# check conservation of mass
plot_conservation_of_mass(
        time_741, mass_from_flux_741, mass_from_con_741,
        "bilayer-741-cons-of-mass.png"
    )

# plot concentration over distance
#Layer thicknesses
W_thickness = 10          # mum
CuEntry_thickness = 65    # mum
Cu_thickness = 175        # mum

# Vector-postprocessor CSV files
W_CSV = glob.glob("val_W_entryCu_Cu_741_W_*.csv")[0]
CUENTRY_CSV = glob.glob("val_W_entryCu_Cu_741_CuEntry_*.csv")[0]
CU_CSV = glob.glob("val_W_entryCu_Cu_741_Cu_*.csv")[0]

X_COL = "x"         
W_VAR_COL = "C_M_W"
CUENTRY_VAR_COL = "C_M_CuEntry"
CU_VAR_COL = "C_M_Cu"

# --- interface x-coordinates (domain order: W -> CuEntry -> Cu) ---
x_W_CuEntry = W_thickness
x_CuEntry_Cu = W_thickness + CuEntry_thickness
x_end = W_thickness + CuEntry_thickness + Cu_thickness

# --- load data directly from files ---
df_W = pd.read_csv(W_CSV)
df_CuEntry = pd.read_csv(CUENTRY_CSV)
df_Cu = pd.read_csv(CU_CSV)

conc_vs_x = {
    "W": (df_W[X_COL].to_numpy(), df_W[W_VAR_COL].to_numpy()),
    "CuEntry": (df_CuEntry[X_COL].to_numpy(), df_CuEntry[CUENTRY_VAR_COL].to_numpy()),
    "Cu": (df_Cu[X_COL].to_numpy(), df_Cu[CU_VAR_COL].to_numpy()),
}

# region shading definitions: (label, x_start, x_end, color)
shading = [
    ("W", 0, x_W_CuEntry, "gray"),
    ("Cu-entry", x_W_CuEntry, x_CuEntry_Cu, "bisque"),
    ("Copper", x_CuEntry_Cu, x_end, "orange"),
]

# --- Combined plot: full system, all three layers ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
for mat in ["W", "CuEntry", "Cu"]:
    x, c = conc_vs_x[mat]
    ax.plot(x, c, c="k")
ax.axvline(x_W_CuEntry, color="b", linestyle="--", linewidth=1.5, label="W\u2013CuEntry interface")
ax.axvline(x_CuEntry_Cu, color="g", linestyle="--", linewidth=1.5, label="CuEntry\u2013Cu interface")
add_shading(ax, shading)
format_ax(ax, (0, x_end))
ax.set_title("Concentration over Distance \u2014 Full System")
plt.savefig("W-CuEntry-Cu-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- W layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
x, c = conc_vs_x["W"]
ax.plot(x, c, c="k")
add_shading(ax, [shading[0]])
format_ax(ax, (0, x_W_CuEntry))
ax.set_title("Concentration over Distance \u2014 W Layer")
plt.savefig("W-CuEntry-Cu-tungsten-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- CuEntry layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
x, c = conc_vs_x["CuEntry"]
ax.plot(x, c, c="k")
add_shading(ax, [shading[1]])
format_ax(ax, (x_W_CuEntry, x_CuEntry_Cu))
ax.set_title("Concentration over Distance \u2014 Cu-Entry Layer")
plt.savefig("W-CuEntry-Cu-entry-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# --- Cu layer only ---
fig, ax = plt.subplots(figsize=(6.5, 5.5))
x, c = conc_vs_x["Cu"]
ax.plot(x, c, c="k")
add_shading(ax, [shading[2]])
format_ax(ax, (x_CuEntry_Cu, x_end))
ax.set_title("Concentration over Distance \u2014 Cu Layer")
plt.savefig("W-CuEntry-Cu-copper-only-vs-dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)


# verify W / CuEntry / Cu system penalty
# Plot solubility ratio and concentration ratio to ensure the jump is properly
# enforced, separately for each interface.
t0 = 100  # choose where concentration stabilizes

for case in cases_bilayer:
    # W / CuEntry interface
    plot_penalty_check(
        case,
        var_key="var_sol_1",
        gold_key="gold_sol_1",
        interface_label="W/CuEntry",
        filename_prefix="bilayer-penalty-W-CuEntry",
    )

    # CuEntry / Cu interface
    plot_penalty_check(
        case,
        var_key="var_sol_2",
        gold_key="gold_sol_2",
        interface_label="CuEntry/Cu",
        filename_prefix="bilayer-penalty-CuEntry-Cu",
    )