import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

num_summation_terms = 1000
N_o = 6.115e28  # lattice density of tungsten
c_o = 0.0001  # dissolved gas atom fraction (-)
D = 2.72e-8 # diffusivity (m^2/s)
l = 2.5e-5 # slab thickness (m)

# Extract data from 'gold' TMAP8 run
# if "/tmap8/doc/" in script_folder.lower():  # if in documentation folder
#     csv_folder = "../../../../test/tests/ver-1dd/gold/ver-1dd_out.csv"
# else:  # if in test folder
#     csv_folder = "./gold/ver-1dd_out.csv"
csv_folder = "./tungsten-ver-1dd_out.csv"
tmap_sol = pd.read_csv(csv_folder)
tmap_time = np.array(tmap_sol["time"])
tmap_prediction = np.array(tmap_sol["scaled_outflux"])
idx = np.where(tmap_time >= 1e-5)[0][0]

# Calculate the breakthrough time from numerical solution
tmap_slope = (tmap_prediction[idx + 1 :] - tmap_prediction[idx:-1]) / (
    tmap_time[idx + 1 :] - tmap_time[idx:-1]
)
tmap_intercept = tmap_time[np.argmax(tmap_slope) + idx] - (
    tmap_prediction[int(np.argmax(tmap_slope) + idx)]
) / np.max(tmap_slope)
output_line = f"The breakthrough time from numerical solution is {tmap_intercept:.4e} s"
print(output_line)

# analytical solution
analytical_time = np.array(tmap_time)
tau_be = l**2 / (2 * (np.pi) ** 2 * D)  # calculate analytical breakthrough time
output_line = f"The breakthrough time from analytical solution is {tau_be:.4e} s"
print(output_line)

def summation_term(num_terms, time):
    sum = 0.0
    for m in range(1, num_terms):
        sum += (-1) ** m * np.exp(-1 * m**2 * time / (2 * tau_be))
    return sum


# Calculate the analytical solution
analytical_flux = (
    N_o * (c_o * D / l) * (1 + 2 * summation_term(num_summation_terms, analytical_time))
)

# Plot figure for verification
x_max = 8*tau_be
y_max = 1.2 * np.max(analytical_flux)
x_tangent_end = 2*tau_be

fig = plt.figure(figsize=[6.5, 5.5])
gs = gridspec.GridSpec(1, 1)
ax = fig.add_subplot(gs[0])

ax.plot(tmap_time, tmap_prediction, label=r"TMAP8", c="tab:gray")  # numerical solution
ax.plot(analytical_time, analytical_flux, label=r"Analytical", c="k", linestyle="--")  # analytical solution
ax.plot(
    [tmap_intercept, x_tangent_end],
    [0, (x_tangent_end - tmap_intercept) * np.max(tmap_slope)],
    label=r"Numerical breakthrough time",
    c="tab:brown",
)
ax.plot(
    [tau_be, tau_be],
    [0,y_max * 0.8],
    label=r"Analytical breakthrough time",
    c="tab:brown",
    linestyle="--",
)
ax.set_xlabel("Time (s)")
ax.set_ylabel("Flux (atom/m$^2$s)")
ax.legend(loc="best")
ax.set_xlim(left=0, right=x_max)
ax.set_ylim(bottom=0, top = y_max)
plt.grid(visible=True, which="major", color="0.65", linestyle="--", alpha=0.3)
RMSE = np.sqrt(np.mean((tmap_prediction - analytical_flux)[idx:] ** 2))
RMSPE = RMSE * 100 / np.mean(analytical_flux[idx:])
ax.text(0.1 * x_max, 0.85 * y_max, "RMSPE = %.3e " % RMSPE + "%", fontweight="bold")
ax.text(0.3* x_max, 0.05 * y_max, "Numerical breakthrough time = %.3e " % tmap_intercept + "s",
    fontweight="bold",
)
ax.text(
    0.3 * x_max, 0.1 * y_max, "Analytical breakthrough time = %.3e " % tau_be + "s",
    fontweight="bold",
)
ax.minorticks_on()
plt.savefig("modified-ver-1dd_comparison_diffusion.png", bbox_inches="tight", dpi=300)
plt.close(fig)
