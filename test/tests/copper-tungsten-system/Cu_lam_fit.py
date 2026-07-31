"""
Fit a polynomial relationship L(T) [Angstrom] for the lattice parameter of
copper as a function of temperature, using data digitized from K. Wang and R. R. Reeber, 
"Thermal expansion of copper," High Temp. Mater. Sci., vol. 35, no. 2, 1996.

NOTE: The T = 1150 K row in the source table (L = 3.64652 A) breaks the
monotonic trend of its neighbors (3.67234 A at 1100 K, 3.68077 A at 1200 K).
By default this script EXCLUDES that point from the fit. Set
EXCLUDE_SUSPECT_POINT = False below if you want to include it as-is.
"""

import numpy as np
import matplotlib.pyplot as plt


# Data: T (K), alpha_L (1e-6 / K, not used in the fit but kept for reference),
# L (Angstrom)
data = [
    (20,    0.19,  3.60300),
    (40,    2.35,  3.60310),
    (60,    5.55,  3.60338),
    (80,    8.31,  3.60388),
    (100,  10.42,  3.60455),
    (150,  13.59,  3.60674),
    (200,  15.19,  3.60935),
    (250,  16.13,  3.61218),
    (300,  16.76,  3.61515),
    (350,  17.14,  3.61821),
    (400,  17.56,  3.62135),
    (450,  17.93,  3.62456),
    (500,  18.27,  3.62784),
    (550,  18.59,  3.63118),
    (600,  18.91,  3.63459),
    (650,  19.22,  3.63806),
    (700,  19.53,  3.64160),
    (750,  19.85,  3.64520),
    (800,  20.17,  3.64887),
    (850,  20.50,  3.65261),
    (900,  20.84,  3.65641),
    (950,  21.18,  3.66029),
    (1000, 21.54,  3.66423),
    (1050, 21.90,  3.66825),
    (1100, 22.28,  3.67234),
    (1150, 22.67,  3.64652),   # <-- suspected typo, see note above
    (1200, 23.08,  3.68077),
    (1250, 23.50,  3.68511),
    (1300, 23.93,  3.68953),
    (1350, 24.39,  3.69405),
]

T_all = np.array([row[0] for row in data], dtype=float)
L_all = np.array([row[2] for row in data], dtype=float)

EXCLUDE_SUSPECT_POINT = True
SUSPECT_T = 1150.0

if EXCLUDE_SUSPECT_POINT:
    mask = T_all != SUSPECT_T
    T_fit = T_all[mask]
    L_fit = L_all[mask]
else:
    T_fit = T_all
    L_fit = L_all


# Fit polynomials of increasing degree and compare goodness of fit
def fit_and_report(degree, T, L):
    coeffs = np.polyfit(T, L, degree)
    poly = np.poly1d(coeffs)
    L_pred = poly(T)
    residuals = L - L_pred
    rmse = np.sqrt(np.mean(residuals**2))
    max_abs_err = np.max(np.abs(residuals))
    return coeffs, poly, rmse, max_abs_err

print(f"{'Degree':>6} | {'RMSE (A)':>12} | {'Max abs err (A)':>16}")
print("-" * 40)
results = {}
for degree in range(1, 6):
    coeffs, poly, rmse, max_err = fit_and_report(degree, T_fit, L_fit)
    results[degree] = (coeffs, poly, rmse, max_err)
    print(f"{degree:>6} | {rmse:>12.6e} | {max_err:>16.6e}")


CHOSEN_DEGREE = 3
coeffs, poly, rmse, max_err = results[CHOSEN_DEGREE]

print(f"\nChosen polynomial degree: {CHOSEN_DEGREE}")
print(f"RMSE: {rmse:.6e} A, Max abs error: {max_err:.6e} A")
print("\nCoefficients (highest power first), L(T) in Angstrom, T in Kelvin:")
for i, c in enumerate(coeffs):
    power = len(coeffs) - 1 - i
    print(f"  T^{power}: {c: .10e}")

# Print as a ready-to-paste Python function
print("\nPython function:\n")
print("def lattice_parameter_Cu(T):")
print('    """Lattice parameter of Cu (Angstrom) as a function of T (K)."""')
terms = " + ".join(
    f"({c:.10e})*T**{len(coeffs)-1-i}" if (len(coeffs)-1-i) != 0 else f"({c:.10e})"
    for i, c in enumerate(coeffs)
)
print(f"    return {terms}")


# Plot data vs fit
T_dense = np.linspace(T_fit.min(), T_fit.max(), 500)
L_dense = poly(T_dense)

plt.figure(figsize=(7, 5))
plt.scatter(T_all, L_all, color="black", s=20, label="Digitized data (incl. suspect pt)")
if EXCLUDE_SUSPECT_POINT:
    plt.scatter([SUSPECT_T], [L_all[T_all == SUSPECT_T][0]],
                color="red", marker="x", s=80, label="Excluded (likely typo)")
plt.plot(T_dense, L_dense, color="tab:blue", label=f"Degree-{CHOSEN_DEGREE} polynomial fit")
plt.xlabel("Temperature (K)")
plt.ylabel("Lattice parameter L (Å)")
plt.title("Cu lattice parameter vs. temperature")
plt.legend()
plt.tight_layout()
plt.savefig("cu_lattice_parameter_fit.png", dpi=150)
plt.show()