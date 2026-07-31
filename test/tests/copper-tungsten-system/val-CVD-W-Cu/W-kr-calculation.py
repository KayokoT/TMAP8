import math
import os
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np

"""
This file computes the surface recombination coefficient from steady-state
permeation data.

Single layer:
    Jp = (D/L) * sqrt(J0/Kr)

Multilayer:
    Jp = S_surface * sqrt(J0/Kr) /
         (L1/(D1*S1) + L2/(D2*S2) + ...)

where
    D = diffusivity
    S = solubility
    L = layer thickness
    J0 = implantation flux
    Kr = surface recombination coefficient

Arrhenius form (standard convention):
    Kr(T) = Kr0 * exp(-Ea / (kB*T))

Fitting strategy
-----------------
This version fits Kr0 and Ea directly to the measured Jp(T) data using scipy.optimize.curve_fit on the
full physical (multilayer) model.

Because Jp spans many orders of magnitude across temperature, the fit is
performed on log(Jp) rather than Jp directly -- this is the standard way to
fit Arrhenius-type data and avoids the fit being dominated entirely by the
highest-temperature (highest-flux) points.
"""

KB_EV = 8.617333262e-5  # eV/K

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(
    SCRIPT_DIR,
    "gold",
    "permflux-vs-temp-experiment.csv",
)

J0 = 3e20   # D/m^2/s

T_COL_CU = "Cu-temp"
JP_COL_CU = "Cu-ss"

T_COL_W = "CVD-W-temp"
JP_COL_W = "CVD-W-ss"

T_COL_WCU = "CVD-W-Cu-temp"
JP_COL_WCU = "CVD-W-Cu-ss"


# diffusivities
def D_Cu_of_T(T):
    return 1.3e-8 * np.exp(-0.27 / (KB_EV * T))


def D_W_of_T(T):
    return 1.9e-7 * np.exp(-0.77 / (KB_EV * T))


# solubilities
def S_Cu_of_T(T):
    S0 = 3.14e24
    Es = 0.57
    return S0 * np.exp(-Es / (KB_EV * T))


def S_W_of_T(T):
    S0 = 1.87e24
    Es = 1.04
    return S0 * np.exp(-Es / (KB_EV * T))


def S_CuEntry_of_T(T):
    S0 = 1.87e24
    Es = 0.6
    return S0 * np.exp(-Es / (KB_EV * T))


# GEOMETRY
W_ONLY_THICKNESS = 100e-6   # standalone W system
CU_ONLY_THICKNESS = 172e-6  # standalone Cu system
# W/Cu system measurements
W_THICKNESS = 10e-6
CU_THICKNESS = 175e-6
CU_ENTRY_THICKNESS = 65e-6


def resistance_of_T(T, layers):
    """Sum of L/(D*S) across layers. Works for scalar or array T."""
    total = 0.0
    for L, D_func, S_func in layers:
        total = total + L / (D_func(T) * S_func(T))
    return total


def Jp_model(T, lnKr0, Ea, layers, surface_solubility):
    """
    Physical forward model (vectorized over T):
        Jp(T) = sqrt(J0 / Kr(T)) / (S1(T) * resistance(T))
    with
        Kr(T) = exp(lnKr0) * exp(-Ea / (kB*T))
        S1 = solubility of the plasma-facing (upstream) layer
    Fitting in terms of lnKr0 (rather than Kr0 directly) keeps curve_fit
    numerically well-behaved, since Kr0 itself can span many orders of
    magnitude depending on material/units.
    """
    Kr = np.exp(lnKr0) * np.exp(-Ea / (KB_EV * T))
    S1 = surface_solubility(T)
    resistance = resistance_of_T(T, layers)
    return np.sqrt(J0 / Kr) / (S1 * resistance)


def log_Jp_model(T, lnKr0, Ea, layers, surface_solubility):
    """log(Jp) version of the model, used for fitting in log-space."""
    return np.log(Jp_model(T, lnKr0, Ea, layers, surface_solubility))


def fit_material(
    label,
    T_col,
    Jp_col,
    layers,
    surface_solubility,
):
    """
    Parameters
    ----------
    layers : list
        [(thickness, diffusivity_function, solubility_function), ...]
    surface_solubility : function
        Solubility of the plasma-facing surface.
    """

    data = pd.read_csv(CSV_PATH)
    data.columns = data.columns.str.strip()
    data = data.dropna(subset=[T_col, Jp_col])

    temps = data[T_col].astype(float).to_numpy()
    Jp_values = data[Jp_col].astype(float).to_numpy()

    if len(temps) < 2:
        raise ValueError(
            f"[{label}] Need at least two valid data points."
        )

    print(f"\n------ {label} ------")

    # --------------------------------------------------------------
    # Pointwise Kr estimate at each T -- diagnostic only, not used in
    # the fit itself. Correct algebra:
    #
    #   Jp = sqrt(J0/Kr) / (S1 * resistance)
    #   =>  Kr = J0 / (Jp * S1 * resistance)^2
    # --------------------------------------------------------------
    for T, Jp in zip(temps, Jp_values):
        resistance = resistance_of_T(T, layers)
        S1 = surface_solubility(T)
        Kr_point = J0 / (Jp * S1 * resistance) ** 2
        print(
            f"T={T:6.1f} K   "
            f"Jp={Jp:.3e}   "
            f"Resistance={resistance:.3e}   "
            f"Kr(point)={Kr_point:.3e}"
        )

    # --------------------------------------------------------------
    # Direct nonlinear least-squares fit of the physical model to
    # log(Jp) data. 
    log_Jp_values = np.log(Jp_values)

    def model_fn(T, lnKr0, Ea):
        return log_Jp_model(T, lnKr0, Ea, layers, surface_solubility)

    # Initial guesses: Ea ~ 1 eV is a typical recombination barrier;
    # back solve lnKr0 from the first data point given that guess.
    Ea_guess = 1.0
    T0 = temps[0]
    resistance0 = resistance_of_T(T0, layers)
    S0 = surface_solubility(T0)
    Kr0_guess = J0 / (Jp_values[0] * S0 * resistance0) ** 2
    lnKr0_guess = math.log(Kr0_guess) - Ea_guess / (KB_EV * T0)
    p0 = [lnKr0_guess, Ea_guess]

    popt, pcov = curve_fit(
        model_fn,
        temps,
        log_Jp_values,
        p0=p0,
        maxfev=20000,
    )
    lnKr0, Ea = popt
    perr = np.sqrt(np.diag(pcov))
    lnKr0_err, Ea_err = perr

    Kr0 = math.exp(lnKr0)
    Kr0_err = Kr0 * lnKr0_err  # propagate error through exp()

    # R^2 computed in log-space, matching the space the fit was done in
    residuals = log_Jp_values - model_fn(temps, lnKr0, Ea)
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((log_Jp_values - np.mean(log_Jp_values)) ** 2)
    r2 = 1.0 if ss_tot == 0 else 1.0 - ss_res / ss_tot

    print(f"\nKr0 = {Kr0:.3e}  (+/- {Kr0_err:.3e})")
    print(f"Ea  = {Ea:.3f} eV  (+/- {Ea_err:.3f} eV)")
    print(f"R^2 (log space) = {r2:.4f}")

    return Kr0, Ea, r2


Kr0_Cu, Ea_Cu, r2_Cu = fit_material(
    "Cu",
    T_COL_CU,
    JP_COL_CU,
    layers=[
        (CU_ONLY_THICKNESS, D_Cu_of_T, S_Cu_of_T),
    ],
    surface_solubility=S_Cu_of_T,
)

Kr0_W, Ea_W, r2_W = fit_material(
    "W",
    T_COL_W,
    JP_COL_W,
    layers=[
        (W_ONLY_THICKNESS, D_W_of_T, S_W_of_T),
    ],
    surface_solubility=S_W_of_T,
)

Kr0_WCu, Ea_WCu, r2_WCu = fit_material(
    "W/Cu",
    T_COL_WCU,
    JP_COL_WCU,
    layers=[
        (W_THICKNESS, D_W_of_T, S_W_of_T),
        (CU_ENTRY_THICKNESS, D_W_of_T, S_CuEntry_of_T),
        (CU_THICKNESS, D_Cu_of_T, S_Cu_of_T),
    ],
    surface_solubility=S_W_of_T,
)