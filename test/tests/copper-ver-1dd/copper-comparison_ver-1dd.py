import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import pandas as pd
import os

# Changes working directory to script directory (for consistent MooseDocs usage)
script_folder = os.path.dirname(__file__)
os.chdir(script_folder)

num_summation_terms = 1000
F = 96485.332 #C/mol, Faraday's constant
T = 1500  # temperature (K)
R = 8.31446261815324 #J/molK, gas constant
lambdaa = 1e-4*(3.602472 + 6.99e-6*T + 1.45e-7*T**2)  # lattice parameter (mum) || lambda is a python keyword
N_o = 4/(lambdaa**3) #lattice density at/mum^3
D_o = 1.74e-6 * 1e12# diffusivity pre-exponential (mum^2/s)
Ed_J_mol = 42000 #diffusion activation energy (J/mol)
c_o = 0.0001  # dissolved gas atom fraction
D = D_o*np.exp(-Ed_J_mol/(R*T))# diffusivity (mum^2/s)
l = 50  # slab thickness (mum)

csv_folder = "./copper-ver-1dd_out.csv"
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
x_max = 20*tau_be
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
plt.savefig("copper-ver-1dd_comparison_diffusion.png", bbox_inches="tight", dpi=300)
plt.close(fig)
