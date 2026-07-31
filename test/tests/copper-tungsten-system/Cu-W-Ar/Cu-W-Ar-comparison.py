import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import glob
import re

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

"""
 This file performs verification checks on each run of Cu-W-Ar.i by plotting the penalty ratio
 and conservation of mass over time. It also creates three separate graphs of the H concentration vs distance  
 for each material, as well as a graph showing all three materials for the designated "run_of_interest."
 Concentration vs time graphs are created for each material for each run.
 """

# === Input Parameters ======
Cu_thickness = 50
W_thickness = 25
Ar_thickness = 50
run_of_interest = "run01" #chosen run to investigate concentration vs distance
t0 = 1

LINESTYLES = ['-', '--', '-.', ':']

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
    simulation_data = pd.read_csv(file_name)
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


def plot_conservation_of_mass_combined(t, datasets, filename, colors):
    """Plot conservation of mass check for multiple systems on one figure.

    Args:
        t (float, ndarray): time array in seconds
        datasets (list of dict): each entry has keys:
            'flux'    (ndarray): accumulated boundary flux in mol
            'mass'    (ndarray): total mass in domain in mol
            't0'      (float):   time before which data is excluded
            'label'   (str):     legend label
        filename (str): output PNG filename
        colors (list of strings): colors to plot, matching length of datasets
    """
    t = np.asarray(t)

    fig, ax = plt.subplots(figsize=[6.5, 5.5])

    for d, c in zip(datasets, colors):
        flux = np.asarray(d['flux'], dtype=float)
        mass = np.asarray(d['mass'], dtype=float)
        t0_local = d['t0']

        idx = np.where(t > t0_local)[0]
        if len(idx) == 0:
            raise ValueError(f"No timesteps found after t = {t0_local} for '{d['label']}'")
        idx = idx[0]

        t_plot = t[idx:]
        percent_diff = np.abs(flux[idx:] / mass[idx:] - 1.0) * 100.0

        ax.plot(t_plot, percent_diff, c=c, label=d['label'])

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Percent difference (%)")
    ax.set_xlim(left=t[0], right=t[-1])
    ax.set_ylim(bottom=0)
    ax.set_title("Conservation of Mass Check")
    ax.legend()
    ax.minorticks_on()
    plt.grid(which="major", color="0.65", linestyle="--", alpha=0.3)
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


def plot_penalty_check(case, run_id, var_key, gold_key, interface_label, filename_prefix):
    """
    Plots concentration-ratio vs solubility-ratio over time for a single interface,
    to verify the penalty method is correctly enforcing the solubility jump there.
    Requires a t0 defined at the beginning of the file, or the point up to which
    the data shouldn't contribute to the RMSPE.

    case: dict with 'time' and the ratio arrays keyed by var_key/gold_key
    run_id: str identifying this run -- used in title/filename
    var_key, gold_key: dict keys into `case` for the concentration-ratio and
                        solubility-ratio time series for this specific interface
    interface_label: e.g. "Cu/W" -- used in title
    filename_prefix: e.g. "Cu-W-penalty"
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
    ax.set_xlim(left = 0)
    ax.set_title(f"Comparison of solubility and concentration ratios ({interface_label}, {run_id})")
    ax.grid(True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.minorticks_on()
    ax.legend(loc="best")

    plt.savefig(f"{filename_prefix}-{run_id}-check.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_concentration_vs_time(material_key, y_max, title, filename, use_linestyle_cycle):
    """Plot H concentration vs time for a single material across all runs.

    Args:
        material_key (str): key into each run's dict in conc_vs_t (e.g. "conc_W")
        y_max (float): maximum concentration across all runs for this material
        title (str): plot title
        filename (str): output PNG filename
        use_linestyle_cycle (bool): whether to cycle line styles per run
            (matches original per-material behavior)
    """
    fig, ax = plt.subplots(figsize=[6.5, 5.5])

    for i, (run, data) in enumerate(conc_vs_t.items()):
        if use_linestyle_cycle:
            ax.plot(data["time"], data[material_key], label=run,
                     linestyle=LINESTYLES[i % len(LINESTYLES)])
        else:
            ax.plot(data["time"], data[material_key], label=run)

    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"Concentration (at/mum$^3$)")
    ax.legend(loc="best")
    ax.set_xlim(left=0, right=t_max)
    ax.set_ylim(bottom=0, top=1.2 * y_max)
    plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
    ax.set_title(title)
    ax.minorticks_on()
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_concentration_vs_distance(materials, xlim, shading_subset, title, filename,
                                    interface_x=None, interface_label=None):
    """Plot H concentration vs distance for one or more material regions.

    Args:
        materials (list of str): material keys (e.g. ["Cu", "W"]) to plot
        xlim (tuple): (left, right) x-axis limits
        shading_subset (list of tuples): region shading entries, see add_shading
        title (str): plot title
        filename (str): output PNG filename
        interface_x (float, optional): x-position of a vertical interface line
        interface_label (str, optional): legend label for the interface line
    """
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    for mat in materials:
        x, c = conc_vs_x[(run_of_interest, mat)]
        mask = regions[mat](x)
        ax.plot(x[mask], c[mask], c="k")
    if interface_x is not None:
        ax.axvline(interface_x, color='g', linestyle='--', linewidth=1.5, label=interface_label)
    add_shading(ax, shading_subset)
    format_ax(ax, xlim)
    ax.set_title(title)
    plt.savefig(filename, bbox_inches="tight", dpi=300)
    plt.close(fig)


#  ====== Input files ======
# pulls main files for plotting concentration over time and sorts them
main_csvs = glob.glob("run[0-9][0-9]_main.csv")

MAIN_COLUMNS = [
    "time", "concentration_at_x_W", "concentration_at_x_Cu", "concentration_at_x_Ar",
    "gold_solubility_ratio", "variable_ratio",
    "Cu_mass_from_concentration", "Cu_mass_from_flux",
    "W_mass_from_concentration", "W_mass_from_flux",
    "Ar_mass_from_concentration", "Ar_mass_from_flux",
    "mass_from_concentration", "mass_from_flux",
]

conc_vs_t = {}
for m in main_csvs:
    run = os.path.splitext(os.path.basename(m))[0]  # run01, run02, etc.
    (time, conc_W, conc_Cu, conc_Ar, gold_sol, var_sol,
     Cu_mass_from_con, Cu_mass_from_flux, W_mass_from_con, W_mass_from_flux,
     Ar_mass_from_con, Ar_mass_from_flux, mass_from_con, mass_from_flux) = read_csv_from_TMAP8(m, MAIN_COLUMNS)

    conc_vs_t[run] = {
        "time": time,
        "conc_W": conc_W,
        "conc_Cu": conc_Cu,
        "conc_Ar": conc_Ar,
        "sol_ratio": gold_sol,
        "conc_ratio": var_sol,
        "Cu_mass_from_con": Cu_mass_from_con,
        "Cu_mass_from_flux": Cu_mass_from_flux,
        "W_mass_from_con": W_mass_from_con,
        "W_mass_from_flux": W_mass_from_flux,
        "Ar_mass_from_con": Ar_mass_from_con,
        "Ar_mass_from_flux": Ar_mass_from_flux,
        "mass_from_con": mass_from_con,
        "mass_from_flux": mass_from_flux,
    }

# files for plotting concentration over distance
files = glob.glob("run*_*_*.csv")

pattern = re.compile(r"(run\d+)_(Cu|W|Ar)_\d+\.csv$")

var_map = {
    "Cu": "C_M_Cu",
    "W": "C_M_W",
    "Ar": "C_M_Ar"
}

# flat structure: (run, mat) -> (mtime, file)
# choose the most recently created file
best = {}

for f in files:
    m = pattern.search(f)
    if not m:
        continue
    run, mat = m.group(1), m.group(2)
    mtime = os.path.getmtime(f)
    key = (run, mat)
    if key not in best or mtime > best[key][0]:
        best[key] = (mtime, f)

conc_vs_x = {}

for (run, mat), (_, f) in best.items():
    var = var_map[mat]
    x, conc = read_csv_from_TMAP8(f, ["x", var])
    conc_vs_x[(run, mat)] = (x, conc)

# determine shared axis limits across all runs
t_max = np.max(np.concatenate([d["time"] for d in conc_vs_t.values()]))
c_max_W = np.max(np.concatenate([d["conc_W"] for d in conc_vs_t.values()]))
c_max_Cu = np.max(np.concatenate([d["conc_Cu"] for d in conc_vs_t.values()]))
c_max_Ar = np.max(np.concatenate([d["conc_Ar"] for d in conc_vs_t.values()]))

# ========= Concentration vs time, one plot per material ===================
plot_concentration_vs_time("conc_W", c_max_W, "Concentration of H in W Layer over Time",
                            "Cu-W-Ar-tungsten-vs-time.png", use_linestyle_cycle=True)
plot_concentration_vs_time("conc_Cu", c_max_Cu, "Concentration of H in Cu layer over Time",
                            "Cu-W-Ar-copper-vs-time.png", use_linestyle_cycle=True)
plot_concentration_vs_time("conc_Ar", c_max_Ar, "Concentration of H in Ar layer over Time",
                            "Cu-W-Ar-argon-vs-time.png", use_linestyle_cycle=False)

# ============ Comparison of concentration as a function of distance ============
# region boundaries
x_CuW = Cu_thickness
x_WAr = Cu_thickness + W_thickness
x_Ar = Cu_thickness + W_thickness + Ar_thickness

# masking functions
regions = {
    "Cu": lambda x: x < x_CuW,
    "W": lambda x: (x >= x_CuW) & (x < x_WAr),
    "Ar": lambda x: x >= x_WAr
}

shading = [
    ("Copper", 0, x_CuW, "forestgreen"),
    ("Tungsten", x_CuW, x_WAr, "palegreen"),
    ("Argon", x_WAr, x_Ar, "salmon")
]

plot_concentration_vs_distance(
    ["Cu", "W"], (0, x_Ar), shading, "Concentration over Distance — Full System",
    "Cu-W-Ar-vs-dist.png", interface_x=x_CuW, interface_label="Cu–W interface"
)
plot_concentration_vs_distance(
    ["Cu"], (0, x_CuW), [shading[0]], "Concentration over Distance — Cu Layer",
    "Cu-W-Ar-copper-only-vs-dist.png"
)
plot_concentration_vs_distance(
    ["W"], (x_CuW, x_WAr), [shading[1]], "Concentration over Distance — W Layer",
    "Cu-W-Ar-tungsten-only-vs-dist.png"
)
plot_concentration_vs_distance(
    ["Ar"], (x_WAr, x_Ar), [shading[2]], "Concentration over Distance — Ar Layer",
    "Cu-W-Ar-argon-only-vs-dist.png"
)

# ============ Verify Penalty ============
# Plot solubility ratio and concentration ratio to ensure the jump is properly enforced
for run_id, case in conc_vs_t.items():
    plot_penalty_check(
        case,
        run_id=run_id.replace("_main", ""),
        var_key="conc_ratio",
        gold_key="sol_ratio",
        interface_label="Cu/W",
        filename_prefix="Cu-W-penalty",
    )

# ==================================== Verify Conservation of Mass ================================================
colors = ["tab:red", "tab:blue", "tab:green", "k"]
for run_id, case in conc_vs_t.items():
    run_label = run_id.replace("_main", "")
    datasets = [
        {'flux': case["Cu_mass_from_flux"], 'mass': case["Cu_mass_from_con"], 't0': 0.2, 'label': 'Cu only'},
        {'flux': case["W_mass_from_flux"], 'mass': case["W_mass_from_con"], 't0': 0.2, 'label': 'W only'},
        {'flux': case["Ar_mass_from_flux"], 'mass': case["Ar_mass_from_con"], 't0': 0.2, 'label': 'Ar only'},
        {'flux': case["mass_from_flux"], 'mass': case["mass_from_con"], 't0': 0.2, 'label': 'Whole system'},
    ]
    plot_conservation_of_mass_combined(
        case["time"], datasets, f"Cu-W-Ar-conservation-of-mass-{run_label}.png", colors
    )