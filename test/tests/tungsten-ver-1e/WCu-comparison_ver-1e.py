import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
from scipy import special
from numpy import sin, cos, tan, sqrt, exp
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

def get_lambdas_analytical(k, l, a):
    # Calculate lambda values for analytical solution
    lambda_range = np.arange(1e-12, 1e2, 1e-5)
    f = 1 / k * sin(lambda_range) * cos(lambda_range * l / a * k)
    g = cos(lambda_range) * sin(lambda_range * l / a * k)
    idx = np.where(np.diff(np.sign(f + g)))
    lambdas = np.expand_dims(lambda_range[idx][::1], axis=0)
    return lambdas


# ========= Comparison of concentration as a function of time in W side ===================

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

csv_folder = "./WCu-ver-1e_csv.csv"
tmap_sol = pd.read_csv(csv_folder)

tmap_time = tmap_sol["time"]
tmap_conc_W = tmap_sol["concentration_at_x_W"]

ax.plot(tmap_time, tmap_conc_W, label=r"TMAP8-W", c="k")

# Analytical parameters
t0 = 0.2 #time offset to avoid singularities
c0 = 50.7079  # concentration at the Cu free surface (moles/m^3)
a = 5e-5  # thickness of the Cu layer (m)
D_Cu = 5.9968958217e-8  # diffusivity in Cu (m^2/s)
D_W = 2.72e-8  # diffusivity in W (m^2/s)
l = 2.5e-5  # thickness of the W layer (m)
x_W = 24e-6  # depth into W layer from Cu/W interface
x_Cu = -1e-6  # depth into Cu layer from Cu/W interface
csv_folder_2 = "./WCu-ver-1e_vector_postproc_line_0142.csv"
k = sqrt(D_Cu / D_W)
# where we compare analytical and numerical model concentration predictions (m)
x1_Cu = x_Cu + a
x2_W = x_W + a
lambdas = get_lambdas_analytical(k, l, a)
t = np.expand_dims(tmap_time, axis=0)

summation = (
    (
        D_Cu * l * sin(lambdas) * sin(k * l / a * lambdas) * (cos(lambdas) - 1)
        + D_W
        * sin(lambdas)
        * (
            k * l * sin(lambdas) * cos(k * l / a * lambdas)
            - a * sin(k * l / a * lambdas)
        )
    )
    / (
        lambdas
        * (a * D_W + l * D_Cu)
        * (np.power(sin(k * l / a * lambdas), 2) + l / a * np.power(sin(lambdas), 2))
    )
    * sin(k * lambdas * (l + a - x2_W) / a)
    * exp(-D_Cu * np.power(lambdas / a, 2) * t.transpose())
)
sums = np.sum(summation, axis=1)

analytical_conc_W = c0 * (D_Cu * (l + a - x2_W) / (l * D_Cu + a * D_W) + 2 * sums)

xmax1 = 0.125
ymax1 = 1.2 * np.max(tmap_conc_W)
idx = np.where(tmap_time >= t0)[0]
RMSE = np.sqrt(np.mean((tmap_conc_W[idx] - analytical_conc_W[idx]) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc_W[idx])
ax.text(0.7*xmax1, 0.05*ymax1, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_time,
    analytical_conc_W,
    label=r"Analytical-W",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)

ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(left = 0, right = xmax1)
ax.set_ylim(bottom = 0, top = ymax1)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
plt.title("Concentration of H in W Layer over Time")
ax.minorticks_on()
plt.savefig("ver-1e_comparison_time.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Closeup of analytical solution ============
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(tmap_time, tmap_conc_W, label=r"TMAP8", c="k")

ax.plot(
    tmap_time,
    analytical_conc_W,
    label=r"Analytical",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)
xmax2 = 0.008
ymax2 = 1.2
ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(0, xmax2)
ax.set_ylim(-0.2, ymax2)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
plt.title("Close-Up of Concentration in W")
ax.minorticks_on()
plt.savefig("ver-1e_comparison_time_closeup.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ========= Comparison of concentration as a function of time in Cu side ===================

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

tmap_conc_Cu = tmap_sol["concentration_at_x_Cu"]

ax.plot(tmap_time, tmap_conc_Cu, label=r"TMAP8-Cu", c="k")


summation = (
    (
        D_Cu * l * np.power(sin(k * l / a * lambdas), 2) * (cos(lambdas) - 1)
        + D_W
        * sin(k * l / a * lambdas)
        * (
            k * l * sin(lambdas) * cos(k * l / a * lambdas)
            - a * sin(k * l / a * lambdas)
        )
    )
    / (
        lambdas
        * (a * D_W + l * D_Cu)
        * (np.power(sin(k * l / a * lambdas), 2) + l / a * np.power(sin(lambdas), 2))
    )
    * sin(lambdas * x1_Cu / a)
    * exp(-D_Cu * np.power(lambdas / a, 2) * t.transpose())
)
sums = np.sum(summation, axis=1)

analytical_conc_Cu = c0 * (
    (D_Cu * l + (a - x1_Cu) * D_W) / (l * D_Cu + a * D_W) + 2 * sums
)
xmax3 = xmax1
ymax3 = 1.2 * np.max(tmap_conc_Cu)
idx = np.where(tmap_time >= t0)[0]
RMSE = np.sqrt(np.mean((tmap_conc_Cu[idx] - analytical_conc_Cu[idx]) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc_Cu[idx])
ax.text(0.7*xmax3, 0.05*ymax3, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_time,
    analytical_conc_Cu,
    label=r"Analytical-Cu",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)

ax.set_xlabel("Time (s)")
ax.set_ylabel(r"Concentration (moles/m$^3$)")
ax.legend(loc="best")
ax.set_xlim(left = 0,right = xmax3)
ax.set_ylim(bottom = 0, top = ymax3)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
plt.title("Concentration of H in Cu layer over Time")
ax.minorticks_on()
plt.savefig("ver-1e_comparison_time_Cu.png", bbox_inches="tight", dpi=300)
plt.close(fig)

# ============ Comparison of concentration as a function of distance ============
fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

tmap_sol = pd.read_csv(csv_folder_2)
tmap_distance = tmap_sol["x"]
tmap_distance_microns = tmap_distance * 1e6
tmap_conc = tmap_sol["u"]
ax.plot(
    tmap_distance_microns,
    tmap_conc,
    label=r"TMAP8",
    c="k",
)

# TMAP 7 Analytical solution
x = tmap_distance
Cu_conc = c0 * ((a - x) * D_W + l * D_Cu) / (l * D_Cu + a * D_W)
W_conc = c0 * ((l + a - x) * D_Cu) / (l * D_Cu + a * D_W)
analytical_conc = (x < a) * Cu_conc + (x >= a) * W_conc
xmax4 = np.max(tmap_distance_microns)
ymax4 = 1.2 *np.max(tmap_conc)
RMSE = np.sqrt(np.mean((tmap_conc - analytical_conc) ** 2))
err_percent = RMSE * 100 / np.mean(analytical_conc)
ax.text(0.1*xmax4, 0.05*ymax4, "RMSPE = %.2f " % err_percent + "%", fontweight="bold")

ax.plot(
    tmap_distance_microns,
    analytical_conc,
    label=r"Analytical",
    c="tab:cyan",
    linestyle="--",
    dashes=(5, 5),
)
ax.set_xlabel("Distance ($\\mu$m)")
ax.set_ylabel(r"Concentration (moles/m$^3$)")
ax.set_xlim(left=0, right = xmax4)
ax.set_ylim(bottom=0, top = ymax4)
ax.legend(loc="best")
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
plt.title("Comparison of Concentration over Distance")
ax.minorticks_on()
plt.savefig("ver-1e_comparison_dist.png", bbox_inches="tight", dpi=300)
plt.close(fig)