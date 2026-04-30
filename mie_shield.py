import sys
from typing import Callable, Iterable, Literal, TypeAlias

import numpy as np
import scipy.integrate
import scipy.stats
import scipy.optimize
import traceback

if not hasattr(scipy.integrate, "trapz"):
    scipy.integrate.trapz = scipy.integrate.trapezoid

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyMieScatt import MieQ
from PySide6.QtCore import Qt, QThread, QLocale, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QGroupBox,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
    QProgressBar,
    QLabel,
    QScrollArea,
    QSpinBox,
    QTabWidget,
)

MATERIALS_DB = [
    ("C", "Углерод (Soot)", 1800.0),
    ("Mg", "Магний (Metal)", 1738.0),
    ("MgCl2", "Хлорид магния", 2320.0),
    ("ZnCl2", "Хлорид цинка", 2907.0),
    ("MgF2", "Фторид магния", 3180.0),
    ("Al4C3", "Карбид алюминия", 2360.0),
    ("Al", "Алюминий (Metal)", 2700.0),
    ("MgO", "Оксид магния", 3580.0),
    ("Al2O3", "Оксид алюминия (корунд)", 3987.0),
    ("CuZn", "Латунь Cu70/Zn30 (Brass)", 8530.0),
]

DEFAULTS = {
    "D_MIN": 0.05,
    "D_MAX": 5.0,
    "D_GEOMETRIC_MEAN": 0.5,
    "SIGMA_GEOMETRIC": 2.0,
    "D_MONO": 1.0,
    "CUSTOM_A": 0.530789,
    "CUSTOM_MU": 1.167607,
    "CUSTOM_SIGMA": 0.373724,
    "WAVELENGTH_MIN": 0.4,
    "WAVELENGTH_MAX": 50.0,
    "WAVELENGTH_STEP": 0.1,
    "POINTS_D": 300,
    "PATH_LENGTH_M": 2.7,
    "INV_D_MIN": 0.01,
    "INV_D_MAX": 50.0,
    "INV_N_SCAN": 1000,
}

ConcMode: TypeAlias = Literal["Массовая", "Числовая"]
DistributionMode: TypeAlias = Literal["monodisperse", "lognormal", "custom"]
InverseWavelengthMode: TypeAlias = Literal["single", "range"]
InverseInputMode: TypeAlias = Literal["mec", "alpha", "tau", "transmittance", "effective_transmittance"]
OptimizationMode: TypeAlias = Literal["window_only", "full"]
OptimizationCriterion: TypeAlias = Literal["mean", "min"]
RIModel: TypeAlias = Callable[[float], complex]
Fractions: TypeAlias = dict[str, float]
DensityMap: TypeAlias = dict[str, float]
ForwardRow: TypeAlias = dict[str, object]

CONC_MASS: ConcMode = "Массовая"
CONC_NUMBER: ConcMode = "Числовая"

DIST_MONODISPERSE: DistributionMode = "monodisperse"
DIST_LOGNORMAL: DistributionMode = "lognormal"
DIST_CUSTOM: DistributionMode = "custom"

INV_WL_SINGLE: InverseWavelengthMode = "single"
INV_WL_RANGE: InverseWavelengthMode = "range"

INV_MEC: InverseInputMode = "mec"
INV_ALPHA: InverseInputMode = "alpha"
INV_TAU: InverseInputMode = "tau"
INV_TRANSMITTANCE: InverseInputMode = "transmittance"
INV_EFFECTIVE_TRANSMITTANCE: InverseInputMode = "effective_transmittance"

OPT_WINDOW_ONLY: OptimizationMode = "window_only"
OPT_FULL: OptimizationMode = "full"
OPT_MEAN: OptimizationCriterion = "mean"
OPT_MIN: OptimizationCriterion = "min"

DENSITY_FALLBACK = 2000.0
TAU_UNDERFLOW_LIMIT = 745.0

def make_wavelengths(min_w, max_w, step):
    if step <= 0 or max_w <= min_w:
        return np.array([])
    num = int(round((max_w - min_w) / step)) + 1
    return np.linspace(min_w, max_w, num)

def _ri_carbon(lam_um: float) -> complex:
    ln_l = np.log(lam_um)
    n = 1.811 + 0.1263 * ln_l + 0.0270 * ln_l**2 + 0.0417 * ln_l**3
    k = 0.5821 + 0.1213 * ln_l + 0.2309 * ln_l**2 - 0.01 * ln_l**3
    return complex(n, k)


def _ri_mg(lam_um: float) -> complex:
    if lam_um < 0.5:
        return complex(0.37, 3.25)
    if lam_um < 1.0:
        n = 0.24 + 0.26 * (1.0 - lam_um) / 0.5
        k = 4.42 + 1.17 * (1.0 - lam_um) / 0.5
        return complex(n, k)
    if lam_um < 5.0:
        return complex(0.10, 8.5 + 4.0 * (lam_um - 1.0) / 4.0)
    if lam_um <= 24.8:
        return complex(0.05, 12.5 + 10.0 * (lam_um - 5.0) / 19.8)
    return complex(0.03, 22.5 + 5.0 * (lam_um - 24.8) / 5.2)


def _ri_mgcl2(lam_um: float) -> complex:
    if lam_um < 2.0:
        return complex(1.675, 0.0)
    if lam_um < 15.0:
        return complex(1.675 - 0.075 * (lam_um - 2) / 13, 0.001 * (lam_um - 2) / 13)
    if lam_um < 18.0:
        return complex(1.55 - 0.3 * (lam_um - 15) / 3, 0.001 + 0.5 * (lam_um - 15) / 3)
    return complex(1.25 - 0.15 * (lam_um - 18) / 12, 0.5 + 0.3 * (lam_um - 18) / 12)


def _ri_zncl2(lam_um: float) -> complex:
    # Anhydrous α-ZnCl2 (tetragonal I-4̄2d, sp.gr. 122). No measured n,k tables exist.
    # Four-oscillator Lorentz model anchored to:
    #   eps_inf = 2.87  (DFPT, Materials Project mp-22909, orientation-averaged)
    #   eps_0   = 5.24  → Σ Δε_j = 2.37 (DFPT, same source; LST satisfied)
    # Phonon frequencies from Angell & Wong, J.Chem.Phys. 53, 2053 (1970)
    # and Janz & James, Spectrochim.Acta A 30, 717 (1974):
    #   305 cm^-1  ν3 antisymm. Zn-Cl stretch (T2)  — main reststrahlen ~33 µm
    #   230 cm^-1  ν1 symm. Zn-Cl stretch (A1)
    #   105 cm^-1  ν2 Cl-Zn-Cl bend
    #    60 cm^-1  network / acoustic-edge modes
    # Cross-check at sodium-D line: √eps_inf = 1.694 vs CRC n_avg = 1.692 (uniaxial+,
    # n_o=1.681, n_e=1.713). Single Lorentz model is used across the full 0.4–50 µm
    # range; in the visible/near-IR all phonon contributions are negligible and
    # n(λ) → √eps_inf with k → 0 to numerical precision.
    nu = 1e4 / lam_um
    eps_inf = 2.87
    osc_nu0 = [305.0, 230.0, 105.0, 60.0]
    osc_de  = [1.60, 0.50, 0.20, 0.07]
    osc_g   = [25.0, 20.0, 30.0, 30.0]
    eps = eps_inf + 0j
    for nu0, de, gj in zip(osc_nu0, osc_de, osc_g):
        eps += de * nu0**2 / (nu0**2 - nu**2 - 1j * nu * gj)
    nc = np.sqrt(eps)
    return complex(abs(nc.real), abs(nc.imag))


def _ri_mgf2(lam_um: float) -> complex:
    if lam_um <= 7.0:
        ls = lam_um**2
        n_sq = (
            1
            + (0.48755108 * ls) / (ls - 0.04338408**2)
            + (0.39875031 * ls) / (ls - 0.09461442**2)
            + (2.3120353 * ls) / (ls - 23.793604**2)
        )
        return complex(np.sqrt(n_sq), 0.0)
    if lam_um < 10.0:
        return complex(1.31 - 0.1 * (lam_um - 7) / 3, 0.05 * (lam_um - 7) / 3)
    if lam_um < 20.0:
        return complex(1.21 - 0.2 * (lam_um - 10) / 10, 0.05 + 0.3 * (lam_um - 10) / 10)
    return complex(1.01 - 0.15 * (lam_um - 20) / 10, 0.35 + 0.4 * (lam_um - 20) / 10)


def _ri_al4c3(lam_um: float) -> complex:
    if lam_um <= 5.0:
        ls = lam_um**2
        n_sq = 1 + (2.8 * ls) / (ls - 0.02**2) + (0.5 * ls) / (ls - 0.1**2) + (1.2 * ls) / (ls - 15.0**2)
        return complex(np.sqrt(n_sq), 0.0)
    if lam_um < 10.0:
        return complex(2.25 - 0.15 * (lam_um - 5) / 5, 0.08 * (lam_um - 5) / 5)
    if lam_um < 20.0:
        return complex(2.10 - 0.25 * (lam_um - 10) / 10, 0.08 + 0.45 * (lam_um - 10) / 10)
    return complex(1.85 - 0.20 * (lam_um - 20) / 10, 0.53 + 0.55 * (lam_um - 20) / 10)


def _ri_al(lam_um: float) -> complex:
    hc_eV_um = 1.23984198
    omega = hc_eV_um / lam_um
    wp = 14.98
    f0, G0 = 0.523, 0.047
    eps = 1.0 - f0 * wp**2 / (omega * (omega + 1j * G0))
    ld_f = [0.227, 0.050, 0.166, 0.030]
    ld_G = [0.333, 0.312, 1.351, 3.382]
    ld_w = [0.162, 1.544, 1.808, 3.473]
    for fj, Gj, wj in zip(ld_f, ld_G, ld_w):
        eps += fj * wp**2 / (wj**2 - omega**2 - 1j * omega * Gj)
    nc = np.sqrt(eps)
    return complex(abs(nc.real), abs(nc.imag))


def _ri_mgo(lam_um: float) -> complex:
    if lam_um <= 5.4:
        ls = lam_um**2
        n_sq = 2.956362 + 0.02195770 / (ls - 0.01428322) - 0.01062387 * ls - 2.04968e-5 * ls**2
        if n_sq < 1.0:
            n_sq = 1.0
        return complex(np.sqrt(n_sq), 0.0)
    nu = 1e4 / lam_um
    eps_inf = 3.014
    osc_nu0 = [384.0, 405.0, 429.0, 590.0]
    osc_f   = [0.20, 1.85, 0.12, 0.10]
    osc_g   = [19.0, 19.0, 19.0, 25.0]
    eps = eps_inf + 0j
    for nu0, fj, gj in zip(osc_nu0, osc_f, osc_g):
        eps += fj * nu0**2 / (nu0**2 - nu**2 - 1j * nu * gj)
    nc = np.sqrt(eps)
    return complex(abs(nc.real), abs(nc.imag))


def _ri_al2o3(lam_um: float) -> complex:
    if lam_um <= 5.0:
        ls = lam_um**2
        B_o = [1.4313493, 0.65054713, 5.3414021]
        C_o = [0.0726631, 0.1193242, 18.028251]
        B_e = [1.5039759, 0.55069141, 6.5927379]
        C_e = [0.0740288, 0.1216529, 20.072248]
        n2_o = 1.0
        n2_e = 1.0
        for Bi, Ci in zip(B_o, C_o):
            n2_o += Bi * ls / (ls - Ci**2)
        for Bi, Ci in zip(B_e, C_e):
            n2_e += Bi * ls / (ls - Ci**2)
        if n2_o < 1.0:
            n2_o = 1.0
        if n2_e < 1.0:
            n2_e = 1.0
        n_avg = np.sqrt((2.0 * n2_o + n2_e) / 3.0)
        return complex(n_avg, 0.0)
    nu = 1e4 / lam_um
    o_einf = 3.064
    o_nu0 = [385.0, 442.0, 569.0, 635.0]
    o_de  = [0.30, 2.70, 3.00, 0.30]
    o_gr  = [0.015, 0.010, 0.020, 0.020]
    eps_o = o_einf + 0j
    for nu0, de, gr in zip(o_nu0, o_de, o_gr):
        gamma = gr * nu0
        eps_o += de * nu0**2 / (nu0**2 - nu**2 - 1j * nu * gamma)
    e_einf = 3.077
    e_nu0 = [400.0, 583.0]
    e_de  = [6.80, 1.70]
    e_gr  = [0.020, 0.035]
    eps_e = e_einf + 0j
    for nu0, de, gr in zip(e_nu0, e_de, e_gr):
        gamma = gr * nu0
        eps_e += de * nu0**2 / (nu0**2 - nu**2 - 1j * nu * gamma)
    eps_poly = (2.0 * eps_o + eps_e) / 3.0
    nc = np.sqrt(eps_poly)
    return complex(abs(nc.real), abs(nc.imag))


def _ri_cuzn(lam_um: float) -> complex:
    cuzn_lam = np.array([
        0.21, 0.22, 0.24, 0.26, 0.28, 0.30, 0.33, 0.36, 0.40,
        0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00,
        1.50, 2.00, 3.00, 5.00, 8.00, 10.0, 20.0, 30.0, 40.0, 50.0,
    ])
    cuzn_n = np.array([
        1.560, 1.490, 1.550, 1.510, 1.460, 1.470, 1.480, 1.450, 1.445,
        1.094, 0.686, 0.527, 0.450, 0.444, 0.446, 0.473, 0.523, 0.603,
        1.044, 1.711, 3.222, 7.097, 13.38, 16.88, 35.82, 54.72, 80.14, 110.5,
    ])
    cuzn_k = np.array([
        1.880, 1.810, 1.750, 1.720, 1.680, 1.700, 1.730, 1.780, 1.805,
        1.829, 2.250, 2.765, 3.253, 3.695, 4.106, 4.890, 5.650, 6.367,
        9.810, 13.10, 19.12, 29.90, 43.34, 51.60, 86.38, 118.7, 148.2, 170.2,
    ])
    n = float(np.interp(lam_um, cuzn_lam, cuzn_n))
    k = float(np.interp(lam_um, cuzn_lam, cuzn_k))
    return complex(n, k)


RI_MODELS: dict[str, RIModel] = {
    "C": _ri_carbon,
    "Mg": _ri_mg,
    "MgCl2": _ri_mgcl2,
    "ZnCl2": _ri_zncl2,
    "MgF2": _ri_mgf2,
    "Al4C3": _ri_al4c3,
    "Al": _ri_al,
    "MgO": _ri_mgo,
    "Al2O3": _ri_al2o3,
    "CuZn": _ri_cuzn,
}


def get_ri(material: str, lam_um: float) -> complex:
    model = RI_MODELS.get(material)
    if model is None:
        return complex(1.5, 0.0)
    return model(lam_um)


def material_density_map() -> DensityMap:
    return {code: rho for (code, _name, rho) in MATERIALS_DB}


def mixture_density(fractions: Fractions, rho_by_code: DensityMap | None = None) -> float:
    densities = material_density_map() if rho_by_code is None else rho_by_code
    return sum(
        float(frac) * float(densities.get(mat_code, DENSITY_FALLBACK))
        for mat_code, frac in fractions.items()
    )


def monodisperse_particle_mass_kg(D_um: float, fractions: Fractions, rho_by_code: DensityMap | None = None) -> float:
    D_m = D_um * 1e-6
    volume_m3 = (np.pi / 6.0) * (D_m ** 3)
    return mixture_density(fractions, rho_by_code) * volume_m3


def distributed_particle_mass_kg(
    diameters_um: np.ndarray,
    pdf_normalized: np.ndarray,
    fractions: Fractions,
    rho_by_code: DensityMap | None = None,
) -> float:
    diameters_m = diameters_um * 1e-6
    volumes_m3 = (np.pi / 6.0) * (diameters_m**3)
    avg_volume_m3 = float(scipy.integrate.trapz(volumes_m3 * pdf_normalized, diameters_um))
    return mixture_density(fractions, rho_by_code) * avg_volume_m3


def resolve_concentration(conc_mode: ConcMode, conc_value: float, avg_mass_kg: float) -> tuple[float, float]:
    if not np.isfinite(conc_value) or conc_value <= 0:
        raise ValueError("Концентрация должна быть > 0.")
    if not np.isfinite(avg_mass_kg) or avg_mass_kg <= 0:
        raise ValueError("Средняя масса частицы некорректна (<=0 или NaN).")

    if conc_mode == CONC_NUMBER:
        num_conc = conc_value
        mass_conc_g = (num_conc * avg_mass_kg) * 1000.0
    else:
        mass_conc_g = conc_value
        num_conc = (mass_conc_g * 1e-3) / avg_mass_kg

    if not np.isfinite(num_conc) or num_conc <= 0:
        raise ValueError("Числовая концентрация некорректна (<=0 или NaN).")
    return num_conc, mass_conc_g


def transmittance_from_tau(tau: float) -> float:
    return 0.0 if tau > TAU_UNDERFLOW_LIMIT else float(np.exp(-tau))


def make_forward_row(
    lam_um: float,
    c_ext_total_mix: float,
    c_ext_parts: dict[str, float],
    num_conc: float,
    mass_conc_g: float,
    path_length_m: float,
) -> ForwardRow:
    c_ext_m2 = c_ext_total_mix * 1e-12
    alpha = num_conc * c_ext_m2
    tau = alpha * path_length_m
    transmittance = transmittance_from_tau(tau)
    mec = alpha / mass_conc_g

    return {
        "wl": float(lam_um),
        "cext_um2": float(c_ext_total_mix),
        "alpha_1m": float(alpha),
        "tau": float(tau),
        "transmittance": float(transmittance),
        "mec_m2g": float(mec),
        "parts": {k: float(v) for k, v in c_ext_parts.items()},
    }


def summarize_forward_rows(rows: list[ForwardRow], materials: Iterable[str]) -> dict[str, object] | None:
    if not rows:
        return None

    mats = list(materials)
    n_rows = len(rows)
    sum_parts = {k: 0.0 for k in mats}

    sum_cext = 0.0
    sum_alpha = 0.0
    sum_tau = 0.0
    sum_transmittance = 0.0
    sum_mec = 0.0

    for row in rows:
        sum_cext += float(row["cext_um2"])
        sum_alpha += float(row["alpha_1m"])
        sum_tau += float(row["tau"])
        sum_transmittance += float(row["transmittance"])
        sum_mec += float(row["mec_m2g"])
        parts = row["parts"]
        for mat_code in mats:
            sum_parts[mat_code] += float(parts[mat_code])

    avg_tau = sum_tau / n_rows
    return {
        "cext_um2": sum_cext / n_rows,
        "alpha_1m": sum_alpha / n_rows,
        "tau": avg_tau,
        "transmittance": sum_transmittance / n_rows,
        "mec_m2g": sum_mec / n_rows,
        "parts": {mat_code: (sum_parts[mat_code] / n_rows) for mat_code in mats},
        "eff_transmittance": transmittance_from_tau(avg_tau),
    }


def inverse_requires_mass_conc(input_mode: InverseInputMode) -> bool:
    return input_mode in (INV_ALPHA, INV_TAU, INV_TRANSMITTANCE, INV_EFFECTIVE_TRANSMITTANCE)


def inverse_uses_transmittance(input_mode: InverseInputMode) -> bool:
    return input_mode in (INV_TRANSMITTANCE, INV_EFFECTIVE_TRANSMITTANCE)


def inverse_metric_label(input_mode: InverseInputMode, use_avg_spectrum: bool) -> str:
    if input_mode == INV_TRANSMITTANCE and use_avg_spectrum:
        return "AVG T"
    if input_mode == INV_EFFECTIVE_TRANSMITTANCE and use_avg_spectrum:
        return "T_eff"
    if inverse_uses_transmittance(input_mode):
        return "T"
    return "AVG MEC" if use_avg_spectrum else "MEC"


def inverse_metric_units(input_mode: InverseInputMode) -> str:
    return "" if inverse_uses_transmittance(input_mode) else "м²/г"


def inverse_solution_names(use_avg_spectrum: bool) -> tuple[str, str]:
    return ("AVG α", "AVG τ") if use_avg_spectrum else ("α", "τ")


def inverse_transmittance_is_reference(input_mode: InverseInputMode, use_avg_spectrum: bool) -> bool:
    return input_mode == INV_TRANSMITTANCE and use_avg_spectrum


def inverse_shows_avg_mec_reference(input_mode: InverseInputMode) -> bool:
    return inverse_uses_transmittance(input_mode)


def resolve_inverse_target(
    input_mode: InverseInputMode,
    target_value: float,
    path_length_m: float,
    mass_conc_g: float | None,
) -> tuple[float, float | None]:
    if inverse_requires_mass_conc(input_mode):
        if mass_conc_g is None or mass_conc_g <= 0:
            raise ValueError("Массовая концентрация должна быть > 0.")
        if inverse_uses_transmittance(input_mode):
            if not (0.0 < target_value <= 1.0):
                raise ValueError("Пропускание T должно быть в диапазоне (0, 1].")
            equivalent_mec = -np.log(target_value) / (mass_conc_g * path_length_m)
            return target_value, equivalent_mec
        if input_mode == INV_TAU:
            return target_value / (mass_conc_g * path_length_m), None
        return target_value / mass_conc_g, None

    return target_value, None


def inverse_metric_from_mec_values(
    mec_values: np.ndarray,
    input_mode: InverseInputMode,
    mass_conc_g: float | None,
    path_length_m: float,
) -> float:
    if input_mode == INV_TRANSMITTANCE:
        if mass_conc_g is None:
            raise ValueError("Массовая концентрация должна быть > 0.")
        tau_vals = np.minimum(mec_values * mass_conc_g * path_length_m, TAU_UNDERFLOW_LIMIT)
        return float(np.mean(np.exp(-tau_vals)))
    if input_mode == INV_EFFECTIVE_TRANSMITTANCE:
        if mass_conc_g is None:
            raise ValueError("Массовая концентрация должна быть > 0.")
        tau_avg = float(np.mean(mec_values) * mass_conc_g * path_length_m)
        return transmittance_from_tau(tau_avg)
    return float(np.mean(mec_values))


def lognormal_pdf(d, d_g, sigma_g):
    d = np.maximum(d, 1e-12)
    return (1.0 / (np.sqrt(2 * np.pi) * d * np.log(sigma_g))) * np.exp(-(np.log(d) - np.log(d_g)) ** 2 / (2 * np.log(sigma_g) ** 2))

def custom_pdf(d, A, mu, sigma):
    d = np.maximum(d, 1e-12)
    return A * (1.0 / (d * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(d) - mu) ** 2 / (2 * sigma ** 2))

def safe_mie_qext(m: complex, lam_nm: float, D_nm: float) -> tuple[float, bool]:
    # Returns (q_ext, ok). q_ext is clamped to >= 0; ok=False on MieQ failure
    # (exception, NaN/inf, or q_ext < -1e-9). Small negative slack q in (-1e-9, 0)
    # is silently clipped to 0 (numerical noise, not a real failure).
    try:
        q_ext = MieQ(m, lam_nm, D_nm, asDict=False)[0]
    except Exception:
        return 0.0, False
    if (not np.isfinite(q_ext)) or (q_ext < -1e-9):
        return 0.0, False
    if q_ext < 0.0:
        q_ext = 0.0
    return float(q_ext), True


def qext_to_cext_um2(q_ext: float, D_nm: float) -> float:
    geom_nm2 = np.pi * (D_nm ** 2) / 4.0
    return q_ext * geom_nm2 / 1e6


def compute_qext_avg(fractions, lam_um, D_um):
    lam_nm = lam_um * 1000.0
    D_nm = D_um * 1000.0
    qext_avg = 0.0
    lost_weight = 0.0
    for mat_code, num_fraction in fractions.items():
        m = get_ri(mat_code, lam_um)
        q_ext, ok = safe_mie_qext(m, lam_nm, D_nm)
        if not ok:
            lost_weight += float(num_fraction)
        qext_avg += float(num_fraction) * q_ext
    return qext_avg, lost_weight

def compute_mec_for_d(D_um, fractions, rho_avg, lam_um):
    D_m = D_um * 1e-6
    qext_avg, lost_weight = compute_qext_avg(fractions, lam_um, D_um)
    if D_m <= 0 or lost_weight > 0.02:
        return np.nan
    mec = (3.0 / 2.0) * (1e-3 / rho_avg) * qext_avg / D_m
    return mec

def nan_safe_bisect(f, a, b, xtol=1e-6, maxiter=100):
    fa = f(a)
    fb = f(b)
    if np.isnan(fa) or np.isnan(fb) or fa * fb >= 0:
        return np.nan
    for _ in range(maxiter):
        if b - a < xtol:
            return (a + b) / 2.0
        mid = (a + b) / 2.0
        fmid = f(mid)
        if np.isnan(fmid):
            left_edge, left_val = a, fa
            for k in range(10, 0, -1):
                probe = a + (mid - a) * k / 11.0
                fp = f(probe)
                if not np.isnan(fp):
                    left_edge, left_val = probe, fp
                    break
            right_edge, right_val = b, fb
            for k in range(1, 11):
                probe = mid + (b - mid) * k / 11.0
                fp = f(probe)
                if not np.isnan(fp):
                    right_edge, right_val = probe, fp
                    break
            if fa * left_val < 0:
                b, fb = left_edge, left_val
            elif right_val * fb < 0:
                a, fa = right_edge, right_val
            elif left_val * right_val < 0:
                a, fa = left_edge, left_val
                b, fb = right_edge, right_val
            else:
                return np.nan
        elif fmid == 0.0:
            return mid
        elif fa * fmid < 0:
            b, fb = mid, fmid
        else:
            a, fa = mid, fmid
    return (a + b) / 2.0


class CalculationWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str)

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.is_aborted = False

    def stop(self):
        self.is_aborted = True

    def run(self):
        p = self.params
        try:
            rho_by_code = material_density_map()
            fractions = p["fractions"]
            is_monodisperse = p.get("monodisperse", False)

            wl_min, wl_max, wl_step = p["wl_range"]
            wavelengths = make_wavelengths(wl_min, wl_max, wl_step)
            if wavelengths.size == 0:
                raise ValueError("Пустая сетка длин волн.")

            actual_step = (wl_max - wl_min) / (wavelengths.size - 1) if wavelengths.size > 1 else wl_step
            self.log_signal.emit(
                f"Спектр: {wavelengths.size} точек, λ=[{wl_min:.4f}–{wl_max:.4f}] мкм, шаг={actual_step:.4f} мкм"
            )
            log_every = max(1, wavelengths.size // 50)
            if log_every > 1:
                self.log_signal.emit(f"ℹ️ Таблица: каждая {log_every}-я точка из {wavelengths.size}")

            conc_mode = p.get("conc_mode", CONC_MASS)
            conc_value = float(p.get("conc_value", 0.01))
            path_length_m = float(p.get("path_length_m", 1.0))
            if not np.isfinite(path_length_m) or path_length_m <= 0:
                raise ValueError("Длина трассы L должна быть > 0.")
            self.log_signal.emit(f"Длина трассы: L={path_length_m:.4f} м")

            if is_monodisperse:
                D_um = float(p["D_um"])

                self.log_signal.emit(f"Монодисперсные частицы: D={D_um:.4f} мкм")

                avg_mass_mixture = monodisperse_particle_mass_kg(D_um, fractions, rho_by_code)
                num_conc, mass_conc_g = resolve_concentration(conc_mode, conc_value, avg_mass_mixture)

                D_nm = D_um * 1000.0

                header = f"{'WL(um)':>10} | {'Cext(um^2)':>15} | {'alpha(1/m)':>15} | {'tau=alpha*L':>15} | {'T':>12} | {'MEC(m^2/g)':>15}"
                self.log_signal.emit(header)
                self.log_signal.emit("-" * len(header))

                results = []

                for i, lam_um in enumerate(wavelengths):
                    if self.is_aborted:
                        self.finished_signal.emit(False, "Остановлено пользователем.")
                        return

                    lam_nm = lam_um * 1000.0
                    c_ext_total_mix = 0.0
                    c_ext_parts = {}
                    mixture_lost_weight = 0.0

                    for mat_code, num_fraction in fractions.items():
                        m = get_ri(mat_code, lam_um)
                        q_ext, ok = safe_mie_qext(m, lam_nm, D_nm)
                        c_ext_um2 = qext_to_cext_um2(q_ext, D_nm)
                        if not ok:
                            mixture_lost_weight += float(num_fraction)

                        c_ext_parts[mat_code] = c_ext_um2
                        c_ext_total_mix += float(num_fraction) * c_ext_um2

                    if mixture_lost_weight > 0.02:
                        msg = f"Критический сбой Mie на {lam_um:.3f} мкм. Потеряно {mixture_lost_weight*100:.1f}% числовой доли смеси."
                        self.finished_signal.emit(False, msg)
                        return

                    row = make_forward_row(lam_um, c_ext_total_mix, c_ext_parts, num_conc, mass_conc_g, path_length_m)
                    results.append(row)

                    if i % log_every == 0 or i == wavelengths.size - 1:
                        self.log_signal.emit(
                            f"{row['wl']:10.4f} | {row['cext_um2']:15.6e} | {row['alpha_1m']:15.6e} | "
                            f"{row['tau']:15.6e} | {row['transmittance']:12.6e} | {row['mec_m2g']:15.6e}"
                        )

                    self.progress_signal.emit(int((i + 1) / wavelengths.size * 100))

                summary = summarize_forward_rows(results, fractions.keys())
                if summary:
                    self.log_signal.emit("-" * len(header))
                    self.log_signal.emit(
                        f"{'AVG':>10} | {summary['cext_um2']:15.6e} | {summary['alpha_1m']:15.6e} | "
                        f"{summary['tau']:15.6e} | {summary['transmittance']:12.6e} | {summary['mec_m2g']:15.6e}"
                    )
                    self.log_signal.emit(f"T_eff=exp(-AVG tau) = {summary['eff_transmittance']:.6e}")

                self.result_signal.emit(
                    {
                        "data": results,
                        "params": p,
                        "num_conc": num_conc,
                        "mass_conc_g": mass_conc_g,
                        "avg_mass_kg": avg_mass_mixture,
                    }
                )
                self.finished_signal.emit(True, "Расчет успешно завершен.")

            else:
                d_min_um = max(p["d_range"][0], 1e-3)
                d_max_um = p["d_range"][1]
                N_D = p.get("points_d", 300)
                dist_type = p.get("dist_type", DIST_LOGNORMAL)

                diameters_um = np.geomspace(d_min_um, d_max_um, N_D)

                if dist_type == DIST_CUSTOM:
                    c_A = p["custom_A"]
                    c_mu = p["custom_mu"]
                    c_sigma = p["custom_sigma"]

                    pdf_values = custom_pdf(diameters_um, c_A, c_mu, c_sigma)
                    mass_numerical = scipy.integrate.trapz(pdf_values, diameters_um)
                    if not np.isfinite(mass_numerical) or mass_numerical <= 1e-12:
                        raise ValueError("Интеграл PDF некорректен (<=0 или NaN).")

                    d_geom_mean = np.exp(c_mu)
                    d_mode_um = np.exp(c_mu - c_sigma**2)

                    self.log_signal.emit(f"Кастомное распределение: A={c_A:.6f}, μ={c_mu:.6f}, σ={c_sigma:.6f}")
                    self.log_signal.emit(f"D_geom_mean=exp(μ)={d_geom_mean:.4f} мкм, Mode={d_mode_um:.4f} мкм")
                    self.log_signal.emit(f"Диапазон: [{d_min_um:.4f}, {d_max_um:.4f}] мкм, N={N_D}")
                    self.log_signal.emit(f"Интеграл (Numerical): {mass_numerical:.6g}")

                    if d_mode_um < d_min_um or d_mode_um > d_max_um:
                        self.log_signal.emit(f"⚠️ ПРЕДУПРЕЖДЕНИЕ: Пик (Mode={d_mode_um:.4f} мкм) вне окна интегрирования!")

                    pdf_normalized = pdf_values / mass_numerical
                    self.log_signal.emit("ℹ️ Модель: Кастомное распределение (перенормировка на 1.0 внутри диапазона).")
                else:
                    dg_um = p["d_dist"][0]
                    sigma_g = p["d_dist"][1]

                    pdf_values = lognormal_pdf(diameters_um, dg_um, sigma_g)
                    mass_numerical = scipy.integrate.trapz(pdf_values, diameters_um)
                    if not np.isfinite(mass_numerical) or mass_numerical <= 1e-12:
                        raise ValueError("Интеграл PDF некорректен (<=0 или NaN).")

                    shape_param = np.log(sigma_g)
                    cdf_max = scipy.stats.lognorm.cdf(d_max_um, shape_param, scale=dg_um)
                    cdf_min = scipy.stats.lognorm.cdf(d_min_um, shape_param, scale=dg_um)
                    coverage_theoretical = cdf_max - cdf_min

                    d_mode_um = dg_um * np.exp(-shape_param**2)

                    self.log_signal.emit(f"Распределение: Dg={dg_um:.4f} мкм, Mode={d_mode_um:.4f} мкм")
                    self.log_signal.emit(f"Диапазон: [{d_min_um:.4f}, {d_max_um:.4f}] мкм, N={N_D}")
                    self.log_signal.emit(f"Покрытие (CDF): {coverage_theoretical*100:.2f}%")
                    self.log_signal.emit(f"Интеграл (Numerical): {mass_numerical:.6g}")

                    if d_mode_um < d_min_um or d_mode_um > d_max_um:
                        self.log_signal.emit(f"⚠️ ПРЕДУПРЕЖДЕНИЕ: Пик (Mode={d_mode_um:.4f} мкм) вне окна интегрирования!")

                    if np.isfinite(coverage_theoretical) and coverage_theoretical < 0.95:
                        self.log_signal.emit("⚠️ ПРЕДУПРЕЖДЕНИЕ: Диапазон отсекает >5% распределения.")

                    if np.isfinite(coverage_theoretical) and coverage_theoretical > 0:
                        rel = abs(mass_numerical - coverage_theoretical) / coverage_theoretical
                        if rel > 0.02:
                            self.log_signal.emit("⚠️ ПРЕДУПРЕЖДЕНИЕ: Численный интеграл отличается от CDF >2%.")

                    pdf_normalized = pdf_values / mass_numerical
                    self.log_signal.emit("ℹ️ Модель: Truncated Log-Normal (перенормировка на 1.0 внутри диапазона).")

                avg_mass_mixture = distributed_particle_mass_kg(diameters_um, pdf_normalized, fractions, rho_by_code)
                num_conc, mass_conc_g = resolve_concentration(conc_mode, conc_value, avg_mass_mixture)

                diameters_nm = diameters_um * 1000.0

                header = f"{'WL(um)':>10} | {'<Cext>(um^2)':>15} | {'alpha(1/m)':>15} | {'tau=alpha*L':>15} | {'T':>12} | {'MEC(m^2/g)':>15}"
                self.log_signal.emit(header)
                self.log_signal.emit("-" * len(header))

                results = []

                for i, lam_um in enumerate(wavelengths):
                    if self.is_aborted:
                        self.finished_signal.emit(False, "Остановлено пользователем.")
                        return

                    lam_nm = lam_um * 1000.0
                    c_ext_total_mix = 0.0
                    c_ext_parts = {}
                    mixture_lost_weight = 0.0

                    for mat_code, num_fraction in fractions.items():
                        m = get_ri(mat_code, lam_um)
                        c_ext_vals_um2 = np.zeros_like(diameters_um)
                        mat_failed_mask = np.zeros_like(diameters_um, dtype=bool)

                        for j, d_nm in enumerate(diameters_nm):
                            q_ext, ok = safe_mie_qext(m, lam_nm, d_nm)
                            if not ok:
                                mat_failed_mask[j] = True
                            c_ext_vals_um2[j] = qext_to_cext_um2(q_ext, d_nm)

                        integrand = c_ext_vals_um2 * pdf_normalized
                        avg_c_ext_mat = float(scipy.integrate.trapz(integrand, diameters_um))

                        c_ext_parts[mat_code] = avg_c_ext_mat
                        c_ext_total_mix += float(num_fraction) * avg_c_ext_mat

                        if np.any(mat_failed_mask):
                            pdf_failed = np.where(mat_failed_mask, pdf_normalized, 0.0)
                            w_lost_mat = float(scipy.integrate.trapz(pdf_failed, diameters_um))
                            mixture_lost_weight += float(num_fraction) * w_lost_mat

                    if mixture_lost_weight > 0.02:
                        msg = f"Критический сбой Mie на {lam_um:.3f} мкм. Потеряно {mixture_lost_weight*100:.1f}% веса распределения."
                        self.finished_signal.emit(False, msg)
                        return

                    row = make_forward_row(lam_um, c_ext_total_mix, c_ext_parts, num_conc, mass_conc_g, path_length_m)
                    results.append(row)

                    if i % log_every == 0 or i == wavelengths.size - 1:
                        self.log_signal.emit(
                            f"{row['wl']:10.4f} | {row['cext_um2']:15.6e} | {row['alpha_1m']:15.6e} | "
                            f"{row['tau']:15.6e} | {row['transmittance']:12.6e} | {row['mec_m2g']:15.6e}"
                        )

                    self.progress_signal.emit(int((i + 1) / wavelengths.size * 100))

                summary = summarize_forward_rows(results, fractions.keys())
                if summary:
                    self.log_signal.emit("-" * len(header))
                    self.log_signal.emit(
                        f"{'AVG':>10} | {summary['cext_um2']:15.6e} | {summary['alpha_1m']:15.6e} | "
                        f"{summary['tau']:15.6e} | {summary['transmittance']:12.6e} | {summary['mec_m2g']:15.6e}"
                    )
                    self.log_signal.emit(f"T_eff=exp(-AVG tau) = {summary['eff_transmittance']:.6e}")

                self.result_signal.emit(
                    {
                        "data": results,
                        "params": p,
                        "num_conc": num_conc,
                        "mass_conc_g": mass_conc_g,
                        "avg_mass_kg": avg_mass_mixture,
                    }
                )
                self.finished_signal.emit(True, "Расчет успешно завершен.")

        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, f"Критическая ошибка:\n{err_msg}")


class InverseWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str)

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.is_aborted = False

    def stop(self):
        self.is_aborted = True

    def run(self):
        p = self.params
        try:
            rho_by_code = material_density_map()
            fractions = p["fractions"]
            rho_avg = mixture_density(fractions, rho_by_code)

            wl_mode = p.get("wl_mode", INV_WL_SINGLE)
            if wl_mode == INV_WL_RANGE:
                wl_min, wl_max, wl_step = p["wl_range"]
                wavelengths = make_wavelengths(wl_min, wl_max, wl_step)
                if wavelengths.size == 0:
                    raise ValueError("Пустая сетка длин волн.")
                actual_wl_step = (wl_max - wl_min) / (wavelengths.size - 1) if wavelengths.size > 1 else wl_step
            else:
                lam_um = float(p["lambda_um"])
                wavelengths = np.array([lam_um], dtype=float)
                actual_wl_step = 0.0
            input_mode = p["input_mode"]
            target_value = float(p["target_value"])
            path_length_m = float(p.get("path_length_m", 1.0))
            if not np.isfinite(path_length_m) or path_length_m <= 0:
                raise ValueError("Длина трассы L должна быть > 0.")

            use_avg_spectrum = wavelengths.size > 1
            metric_label = inverse_metric_label(input_mode, use_avg_spectrum)
            metric_units = inverse_metric_units(input_mode)
            mass_conc_g = float(p["mass_conc_g"]) if inverse_requires_mass_conc(input_mode) else None
            target_mec, equivalent_mec = resolve_inverse_target(input_mode, target_value, path_length_m, mass_conc_g)

            D_min_um = float(p["D_min_um"])
            D_max_um = float(p["D_max_um"])
            N_scan = int(p["N_scan"])

            self.log_signal.emit(f"=== ОБРАТНАЯ ЗАДАЧА (МОНОДИСПЕРС) ===")
            if use_avg_spectrum:
                self.log_signal.emit(
                    f"Спектр: {wavelengths.size} точек, λ=[{wavelengths[0]:.4f}–{wavelengths[-1]:.4f}] мкм, шаг={actual_wl_step:.4f} мкм"
                )
            else:
                self.log_signal.emit(f"λ = {wavelengths[0]:.4f} мкм")
            self.log_signal.emit(f"L = {path_length_m:.4f} м")
            if input_mode == INV_ALPHA:
                alpha_name = "AVG α" if use_avg_spectrum else "α"
                tau_name = "AVG τ" if use_avg_spectrum else "τ"
                self.log_signal.emit(f"{alpha_name} (цель) = {target_value:.6e} 1/м")
                self.log_signal.emit(f"ρ_mass = {mass_conc_g:.6e} г/м³")
                self.log_signal.emit(f"{tau_name} (для этой трассы) = {target_value * path_length_m:.6e}")
            elif input_mode == INV_TAU:
                tau_name = "AVG τ=αL" if use_avg_spectrum else "τ=αL"
                self.log_signal.emit(f"{tau_name} (цель) = {target_value:.6e}")
                self.log_signal.emit(f"ρ_mass = {mass_conc_g:.6e} г/м³")
            elif input_mode == INV_TRANSMITTANCE:
                self.log_signal.emit(f"{metric_label} (цель) = {target_value:.6e}")
                self.log_signal.emit(f"ρ_mass = {mass_conc_g:.6e} г/м³")
                self.log_signal.emit(f"Справочно: -ln(T)/(ρ_mass·L) = {equivalent_mec:.6e} м²/г")
            elif input_mode == INV_EFFECTIVE_TRANSMITTANCE:
                self.log_signal.emit(f"{metric_label}=exp(-AVG τ) (цель) = {target_value:.6e}")
                self.log_signal.emit(f"ρ_mass = {mass_conc_g:.6e} г/м³")
                self.log_signal.emit(f"Эквивалентная AVG MEC = {equivalent_mec:.6e} м²/г")
            units_suffix = f" {metric_units}" if metric_units else ""
            self.log_signal.emit(f"{metric_label} (цель) = {target_mec:.6e}{units_suffix}")
            self.log_signal.emit(f"Диапазон поиска D: [{D_min_um:.4f}, {D_max_um:.4f}] мкм")
            self.log_signal.emit(f"Точек сканирования: {N_scan}")
            self.log_signal.emit(f"ρ_avg (смесь) = {rho_avg:.1f} кг/м³")
            self.log_signal.emit("-" * 60)

            diameters_um = np.geomspace(D_min_um, D_max_um, N_scan)
            mec_values = np.zeros(N_scan)

            def compute_mec_spectrum(D_um):
                vals = np.array([compute_mec_for_d(D_um, fractions, rho_avg, lam) for lam in wavelengths], dtype=float)
                if np.any(~np.isfinite(vals)):
                    return None
                return vals

            def compute_inverse_metric(D_um):
                vals = compute_mec_spectrum(D_um)
                if vals is None:
                    return np.nan
                return inverse_metric_from_mec_values(vals, input_mode, mass_conc_g, path_length_m)

            self.log_signal.emit(f"Сканирование {metric_label}(D)...")

            for i, D_um in enumerate(diameters_um):
                if self.is_aborted:
                    self.finished_signal.emit(False, "Остановлено пользователем.")
                    return

                mec_values[i] = compute_inverse_metric(D_um)

                if i % 100 == 0:
                    self.progress_signal.emit(int(50 * i / N_scan))

            self.progress_signal.emit(50)

            g_values = mec_values - target_mec

            nan_count = np.sum(np.isnan(mec_values))
            if nan_count > 0:
                nan_ranges = []
                in_nan = False
                start_idx = 0
                for i in range(N_scan):
                    if np.isnan(mec_values[i]):
                        if not in_nan:
                            in_nan = True
                            start_idx = i
                    else:
                        if in_nan:
                            nan_ranges.append((diameters_um[start_idx], diameters_um[i-1]))
                            in_nan = False
                if in_nan:
                    nan_ranges.append((diameters_um[start_idx], diameters_um[-1]))

                self.log_signal.emit(f"⚠️ {nan_count} точек с NaN (сбои MieQ) исключены из поиска.")
                for d_start, d_end in nan_ranges[:3]:
                    self.log_signal.emit(f"   NaN диапазон: [{d_start:.4f}, {d_end:.4f}] мкм")
                if len(nan_ranges) > 3:
                    self.log_signal.emit(f"   ... и ещё {len(nan_ranges)-3} диапазонов")

            sign_changes = []
            eps = max(abs(target_mec) * 1e-4, 1e-12)
            eps_near = max(abs(target_mec) * 1e-3, 1e-10)

            interval_indices = set()

            for i in range(N_scan - 1):
                g_i = g_values[i]
                g_i1 = g_values[i + 1]

                if np.isnan(g_i) or np.isnan(g_i1):
                    continue

                if g_i * g_i1 < 0:
                    interval_indices.add(i)

            for i in range(N_scan):
                g_i = g_values[i]

                if np.isnan(g_i):
                    continue

                if abs(g_i) < eps:
                    sign_changes.append((i, "exact"))
                    continue

                if i < N_scan - 1:
                    g_i1 = g_values[i + 1]
                    if not np.isnan(g_i1) and g_i * g_i1 < 0:
                        sign_changes.append((i, "interval"))
                        continue

                if 0 < i < N_scan - 1:
                    if (i in interval_indices) or ((i - 1) in interval_indices):
                        continue

                    g_prev = g_values[i - 1]
                    g_next = g_values[i + 1]
                    if not np.isnan(g_prev) and not np.isnan(g_next):
                        if abs(g_i) < eps_near and abs(g_i) < abs(g_prev) and abs(g_i) < abs(g_next):
                            sign_changes.append((i, "near"))

            self.log_signal.emit(f"Найдено кандидатов на решение: {len(sign_changes)}")

            solutions = []

            for idx, (i, change_type) in enumerate(sign_changes):
                if self.is_aborted:
                    self.finished_signal.emit(False, "Остановлено пользователем.")
                    return

                try:
                    if change_type == "exact" or change_type == "near":
                        D_solution = diameters_um[i]
                    else:
                        D_left = diameters_um[i]
                        D_right = diameters_um[i + 1]

                        if np.isnan(mec_values[i]) or np.isnan(mec_values[i + 1]):
                            continue

                        def objective(D_um):
                            return compute_inverse_metric(D_um) - target_mec

                        D_solution = nan_safe_bisect(objective, D_left, D_right, xtol=1e-6)

                        if np.isnan(D_solution):
                            self.log_signal.emit(f"⚠️ Интервал [{D_left:.4f}, {D_right:.4f}] мкм: bisect не сошёлся (NaN внутри).")
                            continue

                    dedup_threshold = max(1e-4, 1e-3 * D_solution)
                    is_duplicate = False
                    for existing_sol in solutions:
                        if abs(existing_sol["D_um"] - D_solution) < dedup_threshold:
                            is_duplicate = True
                            break

                    if is_duplicate:
                        continue

                    D_m = D_solution * 1e-6
                    V_m3 = (np.pi / 6.0) * (D_m ** 3)
                    m_particle = rho_avg * V_m3

                    mec_check = compute_inverse_metric(D_solution)

                    if np.isnan(mec_check):
                        self.log_signal.emit(f"⚠️ Решение D={D_solution:.4f} мкм отброшено (MieQ сбой).")
                        continue

                    rel_error = abs(mec_check - target_mec) / target_mec if target_mec > 0 else abs(mec_check - target_mec)
                    if rel_error > 0.01:
                        self.log_signal.emit(f"⚠️ Решение D={D_solution:.4f} мкм отброшено (ошибка {rel_error*100:.2f}%).")
                        continue

                    mec_spectrum_check = compute_mec_spectrum(D_solution)
                    avg_mec_check = float(np.mean(mec_spectrum_check)) if mec_spectrum_check is not None else np.nan

                    if inverse_requires_mass_conc(input_mode):
                        N_solution = (mass_conc_g * 1e-3) / m_particle
                        alpha_check = avg_mec_check * mass_conc_g
                        tau_check = alpha_check * path_length_m
                        transmittance_check = transmittance_from_tau(tau_check)
                    else:
                        N_solution = None
                        alpha_check = None
                        tau_check = None
                        transmittance_check = None

                    solutions.append({
                        "D_um": D_solution,
                        "N": N_solution,
                        "mec_check": mec_check,
                        "avg_mec_check": avg_mec_check,
                        "alpha_check": alpha_check,
                        "tau_check": tau_check,
                        "transmittance_check": transmittance_check,
                        "m_particle_kg": m_particle,
                        "rel_error": rel_error,
                    })

                except Exception as e:
                    self.log_signal.emit(f"⚠️ Ошибка на интервале #{i}: {e}")

                self.progress_signal.emit(50 + int(50 * (idx + 1) / max(len(sign_changes), 1)))

            self.log_signal.emit("-" * 60)

            if not solutions:
                self.log_signal.emit("❌ РЕШЕНИЙ НЕ НАЙДЕНО в заданном диапазоне D.")
                self.log_signal.emit("   Попробуйте расширить диапазон или проверить входные данные.")
            else:
                self.log_signal.emit(f"✅ НАЙДЕНО РЕШЕНИЙ: {len(solutions)}")
                self.log_signal.emit("")

                for i, sol in enumerate(solutions):
                    self.log_signal.emit(f"--- Решение #{i+1} ---")
                    self.log_signal.emit(f"  D = {sol['D_um']:.6f} мкм")
                    if sol['N'] is not None:
                        self.log_signal.emit(f"  N = {sol['N']:.6e} 1/м³")
                        alpha_name, tau_name = inverse_solution_names(use_avg_spectrum)
                        transmittance_note = (
                            "справка"
                            if inverse_transmittance_is_reference(input_mode, use_avg_spectrum)
                            else "проверка"
                        )
                        self.log_signal.emit(f"  {alpha_name} (проверка) = {sol['alpha_check']:.6e} 1/м")
                        self.log_signal.emit(f"  {tau_name}=αL (проверка) = {sol['tau_check']:.6e}")
                        self.log_signal.emit(f"  exp(-{tau_name}) ({transmittance_note}) = {sol['transmittance_check']:.6e}")
                    if inverse_shows_avg_mec_reference(input_mode):
                        self.log_signal.emit(f"  AVG MEC (справка) = {sol['avg_mec_check']:.6e} м²/г")
                    self.log_signal.emit(f"  {metric_label} (проверка) = {sol['mec_check']:.6e}{units_suffix}")
                    self.log_signal.emit(f"  Относительная ошибка = {sol['rel_error']*100:.4f}%")
                    self.log_signal.emit(f"  m_particle = {sol['m_particle_kg']:.6e} кг")
                    self.log_signal.emit("")

            self.log_signal.emit("-" * 60)
            self.log_signal.emit(f"Диапазон {metric_label}(D) на сетке:")
            valid_mask = ~np.isnan(mec_values)
            if np.any(valid_mask):
                valid_mec = mec_values[valid_mask]
                valid_D = diameters_um[valid_mask]
                idx_min = np.argmin(valid_mec)
                idx_max = np.argmax(valid_mec)
                self.log_signal.emit(f"  {metric_label}_min = {valid_mec[idx_min]:.6e}{units_suffix} при D = {valid_D[idx_min]:.4f} мкм")
                self.log_signal.emit(f"  {metric_label}_max = {valid_mec[idx_max]:.6e}{units_suffix} при D = {valid_D[idx_max]:.4f} мкм")
            else:
                self.log_signal.emit("  ⚠️ Все точки недействительны (сбои MieQ).")

            self.result_signal.emit({
                "solutions": solutions,
                "scan_D": diameters_um.tolist(),
                "scan_MEC": mec_values.tolist(),
                "target_mec": target_mec,
                "wavelengths": wavelengths.tolist(),
                "metric_label": metric_label,
                "metric_units": metric_units,
                "equivalent_mec": equivalent_mec,
                "params": p,
            })

            self.progress_signal.emit(100)
            self.finished_signal.emit(True, f"Обратная задача завершена. Найдено решений: {len(solutions)}")

        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, f"Критическая ошибка:\n{err_msg}")


class OptimizationWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str)

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.is_aborted = False

    def stop(self):
        self.is_aborted = True

    def run(self):
        p = self.params
        try:
            rho_by_code = material_density_map()
            fractions = p["fractions"]
            rho_avg = mixture_density(fractions, rho_by_code)

            wl_min, wl_max, wl_step = p["wl_range"]
            wavelengths = make_wavelengths(wl_min, wl_max, wl_step)
            n_wl = len(wavelengths)
            if n_wl == 0:
                self.finished_signal.emit(False, "Пустой спектр длин волн.")
                return
            path_length_m = float(p.get("path_length_m", 1.0))
            if not np.isfinite(path_length_m) or path_length_m <= 0:
                raise ValueError("Длина трассы L должна быть > 0.")

            D_scan_min = p["D_scan_min"]
            D_scan_max = p["D_scan_max"]
            N_D_scan = p["N_D_scan"]
            N_D_points = p["N_D_points"]
            N_window_grid = p.get("N_window_grid", 50)
            criterion = p["criterion"]
            mode = p["mode"]
            mu_fixed = p.get("mu_fixed")
            sigma_fixed = p.get("sigma_fixed")
            mu_range = p.get("mu_range")
            sigma_range = p.get("sigma_range")

            D_scan = np.geomspace(D_scan_min, D_scan_max, N_D_scan)

            self.log_signal.emit("=== ОПТИМИЗАЦИЯ MEC ===")
            self.log_signal.emit(f"Спектр: [{wl_min}, {wl_max}] мкм, шаг {wl_step}")
            self.log_signal.emit(f"L = {path_length_m:.4f} м")
            self.log_signal.emit(f"Область D: [{D_scan_min}, {D_scan_max}] мкм, {N_D_scan} точек карты")
            self.log_signal.emit(f"Точек интегрирования PDF: {N_D_points}")
            self.log_signal.emit(f"Критерий: {criterion}")
            self.log_signal.emit(f"Режим: {mode}")
            self.log_signal.emit(f"ρ_avg = {rho_avg:.1f} кг/м³")
            self.log_signal.emit("-" * 60)

            # Phase 1: Build MEC map
            self.log_signal.emit("Фаза 1: Построение карты MEC(D, λ)...")
            mec_map = np.zeros((N_D_scan, n_wl))
            total_phase1 = N_D_scan * n_wl
            done = 0

            for i, D_um in enumerate(D_scan):
                if self.is_aborted:
                    self.finished_signal.emit(False, "Остановлено пользователем.")
                    return
                for j, lam_um in enumerate(wavelengths):
                    mec_map[i, j] = compute_mec_for_d(D_um, fractions, rho_avg, lam_um)
                    done += 1
                pct = int(40 * done / total_phase1)
                self.progress_signal.emit(pct)

            # Send heatmap data
            self.result_signal.emit({
                "type": "heatmap",
                "D_scan": D_scan.tolist(),
                "wavelengths": wavelengths.tolist(),
                "mec_map": mec_map.tolist(),
                "path_score_map": (mec_map * path_length_m).tolist(),
                "path_length_m": path_length_m,
                "metric_label": "MEC·L",
                "metric_units": "м³/г",
            })

            self.log_signal.emit(f"Карта MEC построена: {N_D_scan} x {n_wl}")

            # Build interpolator from log-space D
            from scipy.interpolate import RegularGridInterpolator
            log_D_scan = np.log(D_scan)
            # Replace NaN with 0 for interpolation
            mec_map_clean = np.where(np.isnan(mec_map), 0.0, mec_map)
            interp = RegularGridInterpolator(
                (log_D_scan, wavelengths), mec_map_clean,
                method='linear', bounds_error=False, fill_value=0.0
            )

            # Phase 1.5: Build MEC(D_MIN, D_MAX) landscape
            self.log_signal.emit("Фаза 1.5: Построение карты MEC·L(D_min, D_max)...")
            mu_for_map = mu_fixed if mode == OPT_WINDOW_ONLY else (
                (mu_range[0] + mu_range[1]) / 2.0 if mu_range else 0.0)
            sigma_for_map = sigma_fixed if mode == OPT_WINDOW_ONLY else (
                (sigma_range[0] + sigma_range[1]) / 2.0 if sigma_range else 1.0)

            dmin_grid = np.geomspace(D_scan_min, D_scan_max, N_window_grid)
            dmax_grid = np.geomspace(D_scan_min, D_scan_max, N_window_grid)
            window_mec = np.full((N_window_grid, N_window_grid), np.nan)
            total_window = N_window_grid * N_window_grid
            done_w = 0

            for i, dm in enumerate(dmin_grid):
                if self.is_aborted:
                    self.finished_signal.emit(False, "Остановлено пользователем.")
                    return
                for j, dx in enumerate(dmax_grid):
                    done_w += 1
                    if dx <= dm * 1.01:
                        continue
                    d_g = np.geomspace(dm, dx, N_D_points)
                    pv = custom_pdf(d_g, 1.0, mu_for_map, sigma_for_map)
                    ps = scipy.integrate.trapz(pv, d_g)
                    if ps <= 1e-30:
                        continue
                    pn = pv / ps
                    log_d_g = np.log(d_g)
                    mec_wl = np.zeros(n_wl)
                    for k in range(n_wl):
                        pts = np.column_stack([log_d_g, np.full(N_D_points, wavelengths[k])])
                        mec_wl[k] = scipy.integrate.trapz(interp(pts) * pn, d_g)
                    score_wl = mec_wl * path_length_m
                    if criterion == OPT_MEAN:
                        window_mec[i, j] = np.mean(score_wl)
                    else:
                        window_mec[i, j] = np.min(score_wl)
                pct_w = int(40 + 20 * done_w / total_window)
                self.progress_signal.emit(pct_w)

            self.result_signal.emit({
                "type": "window_heatmap",
                "dmin_grid": dmin_grid.tolist(),
                "dmax_grid": dmax_grid.tolist(),
                "window_mec": window_mec.tolist(),
                "criterion": criterion,
                "path_length_m": path_length_m,
                "metric_label": "MEC·L",
                "metric_units": "м³/г",
            })
            self.log_signal.emit(f"Карта MEC·L(D_min, D_max) построена: {N_window_grid} x {N_window_grid}")

            # Phase 2: Optimization
            self.log_signal.emit("Фаза 2: Оптимизация окна...")

            eval_count = [0]

            def objective(x):
                if self.is_aborted:
                    return 0.0
                if mode == OPT_WINDOW_ONLY:
                    d_min_opt, d_max_opt = x[0], x[1]
                    mu_val = mu_fixed
                    sigma_val = sigma_fixed
                else:
                    d_min_opt, d_max_opt, mu_val, sigma_val = x[0], x[1], x[2], x[3]

                if d_max_opt <= d_min_opt * 1.01:
                    return 1e10

                d_grid = np.geomspace(d_min_opt, d_max_opt, N_D_points)
                pdf_vals = custom_pdf(d_grid, 1.0, mu_val, sigma_val)
                pdf_sum = scipy.integrate.trapz(pdf_vals, d_grid)
                if pdf_sum <= 1e-30:
                    return 1e10
                pdf_norm = pdf_vals / pdf_sum

                log_d_grid = np.log(d_grid)
                mec_per_wl = np.zeros(n_wl)
                for j in range(n_wl):
                    pts = np.column_stack([log_d_grid, np.full(N_D_points, wavelengths[j])])
                    mec_mono = interp(pts)
                    mec_per_wl[j] = scipy.integrate.trapz(mec_mono * pdf_norm, d_grid)

                eval_count[0] += 1
                if eval_count[0] % 50 == 0:
                    pct = min(99, 60 + int(39 * eval_count[0] / 2000))
                    self.progress_signal.emit(pct)

                score_per_wl = mec_per_wl * path_length_m

                if criterion == OPT_MEAN:
                    return -np.mean(score_per_wl)
                else:
                    return -np.min(score_per_wl)

            # Build bounds
            if mode == OPT_WINDOW_ONLY:
                bounds = [
                    (D_scan_min, D_scan_max * 0.9),
                    (D_scan_min * 1.1, D_scan_max),
                ]
            else:
                bounds = [
                    (D_scan_min, D_scan_max * 0.9),
                    (D_scan_min * 1.1, D_scan_max),
                    (mu_range[0], mu_range[1]),
                    (sigma_range[0], sigma_range[1]),
                ]

            result = scipy.optimize.differential_evolution(
                objective, bounds,
                maxiter=200, tol=1e-6, seed=42,
                polish=True, disp=False,
            )

            if self.is_aborted:
                self.finished_signal.emit(False, "Остановлено пользователем.")
                return

            if mode == OPT_WINDOW_ONLY:
                d_min_best, d_max_best = result.x[0], result.x[1]
                mu_best, sigma_best = mu_fixed, sigma_fixed
            else:
                d_min_best, d_max_best = result.x[0], result.x[1]
                mu_best, sigma_best = result.x[2], result.x[3]

            # Compute final MEC spectrum at optimum
            d_grid_best = np.geomspace(d_min_best, d_max_best, N_D_points)
            pdf_best = custom_pdf(d_grid_best, 1.0, mu_best, sigma_best)
            pdf_sum_best = scipy.integrate.trapz(pdf_best, d_grid_best)
            pdf_norm_best = pdf_best / pdf_sum_best if pdf_sum_best > 1e-30 else pdf_best

            mec_spectrum = np.zeros(n_wl)
            log_d_best = np.log(d_grid_best)
            for j in range(n_wl):
                pts = np.column_stack([log_d_best, np.full(N_D_points, wavelengths[j])])
                mec_mono = interp(pts)
                mec_spectrum[j] = scipy.integrate.trapz(mec_mono * pdf_norm_best, d_grid_best)

            mec_mean = float(np.mean(mec_spectrum))
            mec_min = float(np.min(mec_spectrum))
            mec_l_spectrum = mec_spectrum * path_length_m
            mec_l_mean = float(np.mean(mec_l_spectrum))
            mec_l_min = float(np.min(mec_l_spectrum))

            self.log_signal.emit("=" * 60)
            self.log_signal.emit("РЕЗУЛЬТАТ ОПТИМИЗАЦИИ:")
            self.log_signal.emit(f"  D_min* = {d_min_best:.6f} мкм")
            self.log_signal.emit(f"  D_max* = {d_max_best:.6f} мкм")
            self.log_signal.emit(f"  μ*     = {mu_best:.6f}")
            self.log_signal.emit(f"  σ*     = {sigma_best:.6f}")
            self.log_signal.emit(f"  MEC_mean = {mec_mean:.6e} м²/г")
            self.log_signal.emit(f"  MEC_min  = {mec_min:.6e} м²/г")
            self.log_signal.emit(f"  (MEC·L)_mean = {mec_l_mean:.6e} м³/г")
            self.log_signal.emit(f"  (MEC·L)_min  = {mec_l_min:.6e} м³/г")
            self.log_signal.emit(f"  Оценок целевой функции: {eval_count[0]}")
            self.log_signal.emit("=" * 60)

            self.result_signal.emit({
                "type": "result",
                "d_min": d_min_best,
                "d_max": d_max_best,
                "mu": mu_best,
                "sigma": sigma_best,
                "mec_mean": mec_mean,
                "mec_min": mec_min,
                "mec_l_mean": mec_l_mean,
                "mec_l_min": mec_l_min,
                "criterion": criterion,
                "path_length_m": path_length_m,
                "wavelengths": wavelengths.tolist(),
                "mec_spectrum": mec_spectrum.tolist(),
                "mec_l_spectrum": mec_l_spectrum.tolist(),
            })

            self.progress_signal.emit(100)
            self.finished_signal.emit(True, f"Оптимизация завершена. (MEC·L)_mean={mec_l_mean:.4e}, (MEC·L)_min={mec_l_min:.4e}")

        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, f"Критическая ошибка:\n{err_msg}")


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Mie Extinction Calculator (Alpha + MEC)")
        self.resize(1250, 950)
        self.last_results = None
        self.last_inverse_results = None
        self.last_optim_results = None
        self.worker = None
        self.inverse_worker = None
        self.optim_worker = None
        self._heatmap_data = None
        self._window_heatmap_data = None
        self._optim_marker = None

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        left = QVBoxLayout()
        root.addLayout(left, 1)
        right = QVBoxLayout()
        root.addLayout(right, 2)

        self.tabs = QTabWidget()
        left.addWidget(self.tabs)

        tab_forward = QWidget()
        self.tabs.addTab(tab_forward, "Прямая задача")
        self._build_forward_tab(tab_forward)

        tab_inverse = QWidget()
        self.tabs.addTab(tab_inverse, "Обратная задача")
        self._build_inverse_tab(tab_inverse)

        tab_optim = QWidget()
        self.tabs.addTab(tab_optim, "Оптимизация")
        self._build_optim_tab(tab_optim)

        # Matplotlib canvas for optimization heatmaps (in the right panel)
        self.opt_figure = Figure(figsize=(10, 4))
        self.opt_canvas = FigureCanvas(self.opt_figure)
        self.opt_canvas.setMinimumHeight(280)
        self.opt_canvas.setVisible(False)
        right.addWidget(self.opt_canvas)

        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setFont(QFont("Consolas", 9))
        self.log.setLineWrapMode(QTextEdit.NoWrap)
        right.addWidget(self.log)

        self.pbar = QProgressBar()
        right.addWidget(self.pbar)

        self.lbl_st = QLabel("Система готова.")
        right.addWidget(self.lbl_st)

        self.tabs.currentChanged.connect(self._on_tab_changed)

    def _build_forward_tab(self, tab):
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        scroll.setWidget(inner)
        outer.addWidget(scroll)
        layout = QVBoxLayout(inner)

        g1 = QGroupBox("1. Частицы")
        layout.addWidget(g1)
        l1 = QFormLayout(g1)

        self.dist_type_combo = QComboBox()
        self.dist_type_combo.addItem("Монодисперс", DIST_MONODISPERSE)
        self.dist_type_combo.addItem("Лог-нормальное", DIST_LOGNORMAL)
        self.dist_type_combo.addItem("Кастомное", DIST_CUSTOM)
        self.dist_type_combo.setCurrentIndex(self.dist_type_combo.findData(DIST_LOGNORMAL))
        self.dist_type_combo.currentTextChanged.connect(self._on_dist_type_changed)
        l1.addRow("Тип распределения:", self.dist_type_combo)

        self.s_d_mono = self._spin_ndec(DEFAULTS["D_MONO"], 0.001, 100, dec=3, step=0.01)
        self.s_d_mono.setEnabled(False)
        self.lbl_d_mono = QLabel("D (мкм):")
        self.lbl_d_mono.setEnabled(False)
        l1.addRow(self.lbl_d_mono, self.s_d_mono)

        self.s_dmin = self._spin_ndec(DEFAULTS["D_MIN"], 0.001, 100, dec=2, step=0.01)
        self.s_dmax = self._spin_ndec(DEFAULTS["D_MAX"], 0.001, 100, dec=2, step=0.01)
        self.s_dgm = self._spin_ndec(DEFAULTS["D_GEOMETRIC_MEAN"], 0.001, 100, dec=2, step=0.01)
        self.s_sigma = self._spin_ndec(DEFAULTS["SIGMA_GEOMETRIC"], 1.01, 10, dec=2, step=0.01)
        self.s_points = QSpinBox()
        self.s_points.setRange(50, 10000)
        self.s_points.setValue(DEFAULTS["POINTS_D"])

        self.lbl_dmin = QLabel("D min (мкм):")
        self.lbl_dmax = QLabel("D max (мкм):")
        self.lbl_dgm = QLabel("Dg (geom. mean, мкм):")
        self.lbl_sigma = QLabel("Sigma_g (безразм.):")
        self.lbl_points = QLabel("Точек сетки (N):")

        l1.addRow(self.lbl_dmin, self.s_dmin)
        l1.addRow(self.lbl_dmax, self.s_dmax)
        l1.addRow(self.lbl_dgm, self.s_dgm)
        l1.addRow(self.lbl_sigma, self.s_sigma)

        self.s_custom_A = self._spin_ndec(DEFAULTS["CUSTOM_A"], 0.0001, 100, dec=6, step=0.01)
        self.s_custom_mu = self._spin_ndec(DEFAULTS["CUSTOM_MU"], -10, 10, dec=6, step=0.01)
        self.s_custom_sigma = self._spin_ndec(DEFAULTS["CUSTOM_SIGMA"], 0.001, 10, dec=6, step=0.01)
        self.s_custom_A.setEnabled(False)
        self.s_custom_mu.setEnabled(False)
        self.s_custom_sigma.setEnabled(False)

        self.lbl_custom_A = QLabel("A (амплитуда):")
        self.lbl_custom_mu = QLabel("μ (ln-среднее):")
        self.lbl_custom_sigma = QLabel("σ (ln-дисперсия):")
        self.lbl_custom_A.setEnabled(False)
        self.lbl_custom_mu.setEnabled(False)
        self.lbl_custom_sigma.setEnabled(False)

        l1.addRow(self.lbl_custom_A, self.s_custom_A)
        l1.addRow(self.lbl_custom_mu, self.s_custom_mu)
        l1.addRow(self.lbl_custom_sigma, self.s_custom_sigma)

        l1.addRow(self.lbl_points, self.s_points)

        g_wl = QGroupBox("2. Спектр")
        layout.addWidget(g_wl)
        l_wl = QFormLayout(g_wl)
        self.s_wl_min = self._spin_ndec(DEFAULTS["WAVELENGTH_MIN"], 0.1, 1000, dec=2, step=0.01)
        self.s_wl_max = self._spin_ndec(DEFAULTS["WAVELENGTH_MAX"], 0.1, 1000, dec=2, step=0.01)
        self.s_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.1, 100, dec=1, step=0.1)
        l_wl.addRow("Lambda min (мкм):", self.s_wl_min)
        l_wl.addRow("Lambda max (мкм):", self.s_wl_max)
        l_wl.addRow("Шаг (мкм):", self.s_wl_step)

        g2 = QGroupBox("3. Смесь (Численная доля, %)")
        layout.addWidget(g2)
        l2 = QVBoxLayout(g2)
        self.mats = {}
        for c, n, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox(f"{c} - {n}")
            sp = QDoubleSpinBox()
            sp.setLocale(QLocale(QLocale.C))
            sp.setRange(0, 100)
            sp.setDecimals(2)
            sp.setSingleStep(1.0)
            sp.setSuffix("%")
            sp.setAlignment(Qt.AlignRight)
            sp.setEnabled(False)
            if c in ["C", "Mg"]:
                cb.setChecked(True)
                sp.setEnabled(True)
            cb.toggled.connect(sp.setEnabled)
            cb.toggled.connect(lambda ch, s=sp: s.setValue(10) if ch and s.value() == 0 else None)
            h.addWidget(cb, 1)
            h.addWidget(sp, 0)
            l2.addLayout(h)
            self.mats[c] = (cb, sp)

        g_conc = QGroupBox("4. Концентрация и трасса")
        layout.addWidget(g_conc)
        l_conc = QFormLayout(g_conc)
        self.conc_mode = QComboBox()
        self.conc_mode.addItem("Числовая", CONC_NUMBER)
        self.conc_mode.addItem("Массовая", CONC_MASS)
        self.conc_mode.setCurrentIndex(self.conc_mode.findData(CONC_MASS))
        self.conc_value = QDoubleSpinBox()
        self.conc_value.setLocale(QLocale(QLocale.C))
        self.conc_value.setRange(1e-30, 1e18)
        self.conc_value.setDecimals(2)
        self.conc_value.setValue(1.2)
        self.conc_value.setAlignment(Qt.AlignRight)
        self.conc_value.setSuffix(" г/м³")
        self.path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.001, 1000000, dec=3, step=0.1)
        self.path_length.setSuffix(" м")
        l_conc.addRow("Тип:", self.conc_mode)
        l_conc.addRow("Значение:", self.conc_value)
        l_conc.addRow("L трассы:", self.path_length)
        self.conc_mode.currentTextChanged.connect(self._on_conc_mode_changed)

        g3 = QGroupBox("5. Управление")
        layout.addWidget(g3)
        l3 = QVBoxLayout(g3)
        self.btn_run = QPushButton("ЗАПУСК")
        self.btn_run.setFixedHeight(40)
        self.btn_run.clicked.connect(self.start)
        self.btn_stop = QPushButton("СТОП")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop)
        self.btn_save = QPushButton("Экспорт результата")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save)
        l3.addWidget(self.btn_run)
        l3.addWidget(self.btn_stop)
        l3.addWidget(self.btn_save)

        layout.addStretch()

    def _build_inverse_tab(self, tab):
        layout = QVBoxLayout(tab)

        g_input = QGroupBox("1. Входные данные")
        layout.addWidget(g_input)
        l_input = QFormLayout(g_input)

        self.inv_wl_mode = QComboBox()
        self.inv_wl_mode.addItem("Одна λ", INV_WL_SINGLE)
        self.inv_wl_mode.addItem("Диапазон λ", INV_WL_RANGE)
        self.inv_wl_mode.currentTextChanged.connect(self._on_inv_wl_mode_changed)
        l_input.addRow("Спектр:", self.inv_wl_mode)

        self.inv_lambda = self._spin_ndec(0.59, 0.1, 1000, dec=2, step=0.01)
        self.lbl_inv_lambda = QLabel("λ (мкм):")
        l_input.addRow(self.lbl_inv_lambda, self.inv_lambda)

        self.inv_wl_min = self._spin_ndec(0.4, 0.1, 1000, dec=2, step=0.01)
        self.inv_wl_max = self._spin_ndec(0.78, 0.1, 1000, dec=2, step=0.01)
        self.inv_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.001, 100, dec=3, step=0.01)
        self.lbl_inv_wl_min = QLabel("λ min (мкм):")
        self.lbl_inv_wl_max = QLabel("λ max (мкм):")
        self.lbl_inv_wl_step = QLabel("Шаг λ (мкм):")
        l_input.addRow(self.lbl_inv_wl_min, self.inv_wl_min)
        l_input.addRow(self.lbl_inv_wl_max, self.inv_wl_max)
        l_input.addRow(self.lbl_inv_wl_step, self.inv_wl_step)
        self._on_inv_wl_mode_changed(self.inv_wl_mode.currentText())

        self.inv_path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.001, 1000000, dec=3, step=0.1)
        self.inv_path_length.setSuffix(" м")
        l_input.addRow("L трассы:", self.inv_path_length)

        self.inv_input_mode = QComboBox()
        self.inv_input_mode.addItem("T_eff = exp(-AVG alpha*L)", INV_EFFECTIVE_TRANSMITTANCE)
        self.inv_input_mode.addItem("AVG T = mean(exp(-alpha*L))", INV_TRANSMITTANCE)
        self.inv_input_mode.addItem("tau = alpha*L", INV_TAU)
        self.inv_input_mode.addItem("alpha (1/м)", INV_ALPHA)
        self.inv_input_mode.addItem("MEC (м²/г)", INV_MEC)
        self.inv_input_mode.currentTextChanged.connect(self._on_inv_mode_changed)
        l_input.addRow("Тип входа:", self.inv_input_mode)

        self.inv_target_value = QDoubleSpinBox()
        self.inv_target_value.setLocale(QLocale(QLocale.C))
        self.inv_target_value.setRange(1e-30, 1e30)
        self.inv_target_value.setDecimals(6)
        self.inv_target_value.setValue(1e-3)
        self.inv_target_value.setAlignment(Qt.AlignRight)
        l_input.addRow("Значение:", self.inv_target_value)

        self.inv_mass_conc = self._spin_ndec(1.2, 1e-30, 1e18, dec=6, step=0.001)
        self.inv_mass_conc.setSuffix(" г/м³")
        self.lbl_inv_mass_conc = QLabel("ρ_mass (г/м³):")
        l_input.addRow(self.lbl_inv_mass_conc, self.inv_mass_conc)
        self._on_inv_mode_changed(self.inv_input_mode.currentText())

        g_search = QGroupBox("2. Диапазон поиска D")
        layout.addWidget(g_search)
        l_search = QFormLayout(g_search)

        self.inv_d_min = self._spin_ndec(DEFAULTS["INV_D_MIN"], 0.001, 1000, dec=3, step=0.01)
        self.inv_d_max = self._spin_ndec(DEFAULTS["INV_D_MAX"], 0.001, 1000, dec=3, step=1.0)
        self.inv_n_scan = QSpinBox()
        self.inv_n_scan.setRange(100, 100000)
        self.inv_n_scan.setValue(DEFAULTS["INV_N_SCAN"])

        l_search.addRow("D min (мкм):", self.inv_d_min)
        l_search.addRow("D max (мкм):", self.inv_d_max)
        l_search.addRow("Точек сканирования:", self.inv_n_scan)

        g_mat = QGroupBox("3. Смесь (Численная доля, %)")
        layout.addWidget(g_mat)
        l_mat = QVBoxLayout(g_mat)
        self.inv_mats = {}
        for c, n, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox(f"{c} - {n}")
            sp = QDoubleSpinBox()
            sp.setLocale(QLocale(QLocale.C))
            sp.setRange(0, 100)
            sp.setDecimals(2)
            sp.setSingleStep(1.0)
            sp.setSuffix("%")
            sp.setAlignment(Qt.AlignRight)
            sp.setEnabled(False)
            if c == "C":
                cb.setChecked(True)
                sp.setEnabled(True)
                sp.setValue(100)
            cb.toggled.connect(sp.setEnabled)
            cb.toggled.connect(lambda ch, s=sp: s.setValue(10) if ch and s.value() == 0 else None)
            h.addWidget(cb, 1)
            h.addWidget(sp, 0)
            l_mat.addLayout(h)
            self.inv_mats[c] = (cb, sp)

        g_ctrl = QGroupBox("4. Управление")
        layout.addWidget(g_ctrl)
        l_ctrl = QVBoxLayout(g_ctrl)

        self.btn_inv_run = QPushButton("ПОИСК D")
        self.btn_inv_run.setFixedHeight(40)
        self.btn_inv_run.clicked.connect(self.start_inverse)
        self.btn_inv_stop = QPushButton("СТОП")
        self.btn_inv_stop.setEnabled(False)
        self.btn_inv_stop.clicked.connect(self.stop_inverse)
        self.btn_inv_save = QPushButton("Экспорт результата")
        self.btn_inv_save.setEnabled(False)
        self.btn_inv_save.clicked.connect(self.save_inverse)

        l_ctrl.addWidget(self.btn_inv_run)
        l_ctrl.addWidget(self.btn_inv_stop)
        l_ctrl.addWidget(self.btn_inv_save)

        layout.addStretch()

    def _build_optim_tab(self, tab):
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        inner = QWidget()
        scroll.setWidget(inner)
        outer.addWidget(scroll)
        layout = QVBoxLayout(inner)

        g_wl = QGroupBox("1. Спектр, трасса и материалы")
        layout.addWidget(g_wl)
        l_wl = QFormLayout(g_wl)
        self.opt_wl_min = self._spin_ndec(DEFAULTS["WAVELENGTH_MIN"], 0.1, 1000, dec=2, step=0.01)
        self.opt_wl_max = self._spin_ndec(DEFAULTS["WAVELENGTH_MAX"], 0.1, 1000, dec=2, step=0.01)
        self.opt_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.1, 100, dec=1, step=0.1)
        l_wl.addRow("λ min (мкм):", self.opt_wl_min)
        l_wl.addRow("λ max (мкм):", self.opt_wl_max)
        l_wl.addRow("Шаг λ (мкм):", self.opt_wl_step)
        self.opt_path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.001, 1000000, dec=3, step=0.1)
        self.opt_path_length.setSuffix(" м")
        l_wl.addRow("L трассы:", self.opt_path_length)

        g_mat = QGroupBox("2. Смесь (Численная доля, %)")
        layout.addWidget(g_mat)
        l_mat = QVBoxLayout(g_mat)
        self.opt_mats = {}
        for c, n, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox(f"{c} - {n}")
            sp = QDoubleSpinBox()
            sp.setLocale(QLocale(QLocale.C))
            sp.setRange(0, 100)
            sp.setDecimals(2)
            sp.setSingleStep(1.0)
            sp.setSuffix("%")
            sp.setAlignment(Qt.AlignRight)
            sp.setEnabled(False)
            if c == "C":
                cb.setChecked(True)
                sp.setEnabled(True)
                sp.setValue(100)
            cb.toggled.connect(sp.setEnabled)
            cb.toggled.connect(lambda ch, s=sp: s.setValue(10) if ch and s.value() == 0 else None)
            h.addWidget(cb, 1)
            h.addWidget(sp, 0)
            l_mat.addLayout(h)
            self.opt_mats[c] = (cb, sp)

        g_dist = QGroupBox("3. Параметры распределения")
        layout.addWidget(g_dist)
        l_dist = QFormLayout(g_dist)
        self.opt_mode = QComboBox()
        self.opt_mode.addItem("Только окно D", OPT_WINDOW_ONLY)
        self.opt_mode.addItem("Полная (μ, σ, окно)", OPT_FULL)
        self.opt_mode.currentTextChanged.connect(self._on_optim_mode_changed)
        l_dist.addRow("Режим:", self.opt_mode)

        self.opt_mu = self._spin_ndec(DEFAULTS["CUSTOM_MU"], -10, 10, dec=6, step=0.01)
        self.opt_sigma = self._spin_ndec(DEFAULTS["CUSTOM_SIGMA"], 0.001, 10, dec=6, step=0.01)
        self.lbl_opt_mu = QLabel("μ (фикс.):")
        self.lbl_opt_sigma = QLabel("σ (фикс.):")
        l_dist.addRow(self.lbl_opt_mu, self.opt_mu)
        l_dist.addRow(self.lbl_opt_sigma, self.opt_sigma)

        self.opt_mu_min = self._spin_ndec(-2.0, -10, 10, dec=6, step=0.1)
        self.opt_mu_max = self._spin_ndec(3.0, -10, 10, dec=6, step=0.1)
        self.opt_sigma_min = self._spin_ndec(0.1, 0.001, 10, dec=6, step=0.01)
        self.opt_sigma_max = self._spin_ndec(2.0, 0.001, 10, dec=6, step=0.01)
        self.lbl_opt_mu_min = QLabel("μ min:")
        self.lbl_opt_mu_max = QLabel("μ max:")
        self.lbl_opt_sigma_min = QLabel("σ min:")
        self.lbl_opt_sigma_max = QLabel("σ max:")
        l_dist.addRow(self.lbl_opt_mu_min, self.opt_mu_min)
        l_dist.addRow(self.lbl_opt_mu_max, self.opt_mu_max)
        l_dist.addRow(self.lbl_opt_sigma_min, self.opt_sigma_min)
        l_dist.addRow(self.lbl_opt_sigma_max, self.opt_sigma_max)
        for w in [self.opt_mu_min, self.opt_mu_max, self.opt_sigma_min, self.opt_sigma_max,
                  self.lbl_opt_mu_min, self.lbl_opt_mu_max, self.lbl_opt_sigma_min, self.lbl_opt_sigma_max]:
            w.setEnabled(False)

        g_search = QGroupBox("4. Область поиска D")
        layout.addWidget(g_search)
        l_search = QFormLayout(g_search)
        self.opt_d_scan_min = self._spin_ndec(0.01, 0.001, 100, dec=3, step=0.01)
        self.opt_d_scan_max = self._spin_ndec(20.0, 0.01, 1000, dec=3, step=1.0)
        self.opt_n_d_scan = QSpinBox()
        self.opt_n_d_scan.setRange(20, 1000)
        self.opt_n_d_scan.setValue(150)
        self.opt_n_d_points = QSpinBox()
        self.opt_n_d_points.setRange(20, 2000)
        self.opt_n_d_points.setValue(200)
        self.opt_n_window = QSpinBox()
        self.opt_n_window.setRange(10, 200)
        self.opt_n_window.setValue(50)
        l_search.addRow("D scan min (мкм):", self.opt_d_scan_min)
        l_search.addRow("D scan max (мкм):", self.opt_d_scan_max)
        l_search.addRow("Точек карты (N_D_scan):", self.opt_n_d_scan)
        l_search.addRow("Точек PDF (N_D_points):", self.opt_n_d_points)
        l_search.addRow("Сетка окна (N_window):", self.opt_n_window)

        g_crit = QGroupBox("5. Критерий")
        layout.addWidget(g_crit)
        l_crit = QFormLayout(g_crit)
        self.opt_criterion = QComboBox()
        self.opt_criterion.addItem("Максимизировать mean MEC·L", OPT_MEAN)
        self.opt_criterion.addItem("Максимизировать min MEC·L", OPT_MIN)
        l_crit.addRow("Критерий:", self.opt_criterion)

        g_ctrl = QGroupBox("6. Управление")
        layout.addWidget(g_ctrl)
        l_ctrl = QVBoxLayout(g_ctrl)
        self.btn_opt_run = QPushButton("ОПТИМИЗАЦИЯ")
        self.btn_opt_run.setFixedHeight(40)
        self.btn_opt_run.clicked.connect(self.start_optim)
        self.btn_opt_stop = QPushButton("СТОП")
        self.btn_opt_stop.setEnabled(False)
        self.btn_opt_stop.clicked.connect(self.stop_optim)
        self.btn_opt_save = QPushButton("Экспорт результата")
        self.btn_opt_save.setEnabled(False)
        self.btn_opt_save.clicked.connect(self.save_optim)
        l_ctrl.addWidget(self.btn_opt_run)
        l_ctrl.addWidget(self.btn_opt_stop)
        l_ctrl.addWidget(self.btn_opt_save)

        layout.addStretch()

    def _on_optim_mode_changed(self, _text=None):
        opt_mode = self._combo_data(self.opt_mode, OPT_WINDOW_ONLY)
        is_full = (opt_mode == OPT_FULL)
        self.opt_mu.setEnabled(not is_full)
        self.opt_sigma.setEnabled(not is_full)
        self.lbl_opt_mu.setEnabled(not is_full)
        self.lbl_opt_sigma.setEnabled(not is_full)
        for w in [self.opt_mu_min, self.opt_mu_max, self.opt_sigma_min, self.opt_sigma_max,
                  self.lbl_opt_mu_min, self.lbl_opt_mu_max, self.lbl_opt_sigma_min, self.lbl_opt_sigma_max]:
            w.setEnabled(is_full)

    def _on_tab_changed(self, index):
        is_optim = (self.tabs.tabText(index) == "Оптимизация")
        self.opt_canvas.setVisible(is_optim)

    def start_optim(self):
        if self.optim_worker and self.optim_worker.isRunning():
            return

        fracs = {}
        tot = 0.0
        for c, (cb, sp) in self.opt_mats.items():
            if cb.isChecked():
                fracs[c] = float(sp.value())
                tot += float(sp.value())
        if not fracs or tot <= 0:
            QMessageBox.warning(self, "Ошибка", "Не выбран состав смеси.")
            return
        norm_fracs = {k: v / tot for k, v in fracs.items()}

        wl_min = float(self.opt_wl_min.value())
        wl_max = float(self.opt_wl_max.value())
        wl_step = float(self.opt_wl_step.value())
        if wl_min >= wl_max or wl_step <= 0:
            QMessageBox.warning(self, "Ошибка", "Неверные параметры спектра.")
            return
        path_length_m = float(self.opt_path_length.value())
        if path_length_m <= 0 or (not np.isfinite(path_length_m)):
            QMessageBox.warning(self, "Ошибка", "Длина трассы L должна быть > 0.")
            return

        d_scan_min = float(self.opt_d_scan_min.value())
        d_scan_max = float(self.opt_d_scan_max.value())
        if d_scan_min >= d_scan_max:
            QMessageBox.warning(self, "Ошибка", "D_scan_min >= D_scan_max.")
            return

        opt_mode = self._combo_data(self.opt_mode, OPT_WINDOW_ONLY)
        is_full = (opt_mode == OPT_FULL)
        criterion = self._combo_data(self.opt_criterion, OPT_MEAN)

        p = {
            "fractions": norm_fracs,
            "wl_range": (wl_min, wl_max, wl_step),
            "path_length_m": path_length_m,
            "D_scan_min": d_scan_min,
            "D_scan_max": d_scan_max,
            "N_D_scan": int(self.opt_n_d_scan.value()),
            "N_D_points": int(self.opt_n_d_points.value()),
            "N_window_grid": int(self.opt_n_window.value()),
            "criterion": criterion,
            "mode": OPT_FULL if is_full else OPT_WINDOW_ONLY,
        }

        if is_full:
            p["mu_range"] = (float(self.opt_mu_min.value()), float(self.opt_mu_max.value()))
            p["sigma_range"] = (float(self.opt_sigma_min.value()), float(self.opt_sigma_max.value()))
            p["mu_fixed"] = None
            p["sigma_fixed"] = None
        else:
            p["mu_fixed"] = float(self.opt_mu.value())
            p["sigma_fixed"] = float(self.opt_sigma.value())

        self.log.clear()
        self.optim_worker = OptimizationWorker(p)
        self.optim_worker.log_signal.connect(self.log_msg)
        self.optim_worker.progress_signal.connect(self.pbar.setValue)
        self.optim_worker.result_signal.connect(self.on_optim_result)
        self.optim_worker.finished_signal.connect(self.on_optim_finish)

        self.btn_opt_run.setEnabled(False)
        self.btn_opt_stop.setEnabled(True)
        self.btn_opt_save.setEnabled(False)
        self.lbl_st.setText("Оптимизация...")
        self.pbar.setValue(0)
        self.optim_worker.start()

    def stop_optim(self):
        if self.optim_worker:
            self.optim_worker.stop()
            self.lbl_st.setText("Прерывание...")

    def on_optim_result(self, res):
        rtype = res.get("type")
        if rtype == "heatmap":
            self._heatmap_data = res
            self._window_heatmap_data = None
            self._optim_marker = None
            self._draw_optim_plots()
        elif rtype == "window_heatmap":
            self._window_heatmap_data = res
            self._draw_optim_plots()
        elif rtype == "result":
            self.last_optim_results = res
            self.btn_opt_save.setEnabled(True)
            marker_value = (
                res.get("mec_l_mean", res["mec_mean"])
                if res.get("criterion") == OPT_MEAN
                else res.get("mec_l_min", res["mec_min"])
            )
            self._optim_marker = (res["d_min"], res["d_max"],
                                  marker_value)
            self._draw_optim_plots()

    def _draw_optim_plots(self):
        self.opt_figure.clear()
        # Disconnect previous motion handler if any
        if hasattr(self, '_opt_motion_cid') and self._opt_motion_cid is not None:
            self.opt_canvas.mpl_disconnect(self._opt_motion_cid)
            self._opt_motion_cid = None

        has_left = getattr(self, '_heatmap_data', None) is not None
        has_right = getattr(self, '_window_heatmap_data', None) is not None
        ncols = (1 if has_left else 0) + (1 if has_right else 0)
        if ncols == 0:
            return
        col = 1

        self._opt_ax_left = None
        self._opt_ax_right = None

        if has_left:
            res = self._heatmap_data
            D_scan = np.array(res["D_scan"])
            wavelengths = np.array(res["wavelengths"])
            metric_map = np.array(res.get("path_score_map", res["mec_map"]))
            metric_label = res.get("metric_label", "MEC")
            metric_units = res.get("metric_units", "м²/г")
            ax1 = self.opt_figure.add_subplot(1, ncols, col)
            col += 1
            mec_plot = np.where(metric_map > 0, metric_map, np.nan)
            im1 = ax1.pcolormesh(wavelengths, D_scan, mec_plot, shading='auto', cmap='hot')
            ax1.set_yscale('log')
            ax1.set_xlabel('λ (мкм)')
            ax1.set_ylabel('D (мкм)')
            ax1.set_title(f'{metric_label}(D, λ) {metric_units}')
            self.opt_figure.colorbar(im1, ax=ax1)
            self._opt_ax_left = ax1
            self._opt_left_data = (wavelengths, D_scan, metric_map, metric_label, metric_units)

        if has_right:
            res_w = self._window_heatmap_data
            dmin_grid = np.array(res_w["dmin_grid"])
            dmax_grid = np.array(res_w["dmax_grid"])
            window_mec = np.array(res_w["window_mec"])
            metric_label = res_w.get("metric_label", "MEC")
            metric_units = res_w.get("metric_units", "м²/г")
            crit_label = "mean" if res_w["criterion"] == OPT_MEAN else "min"

            ax2 = self.opt_figure.add_subplot(1, ncols, col)
            im2 = ax2.pcolormesh(dmax_grid, dmin_grid, window_mec, shading='auto', cmap='viridis')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.set_xlabel('D_max (мкм)')
            ax2.set_ylabel('D_min (мкм)')
            ax2.set_title(f'{metric_label}_{crit_label}(D_min, D_max)')
            self.opt_figure.colorbar(im2, ax=ax2)
            self._opt_ax_right = ax2
            self._opt_right_data = (dmax_grid, dmin_grid, window_mec, metric_label, metric_units)

            marker = getattr(self, '_optim_marker', None)
            if marker is not None:
                d_min_best, d_max_best, mec_best = marker
                ax2.plot(d_max_best, d_min_best, marker='*', color='red',
                         markersize=15, markeredgecolor='white', markeredgewidth=1.0)
                ax2.annotate(f'{mec_best:.3e}',
                             xy=(d_max_best, d_min_best),
                             xytext=(8, 8), textcoords='offset points',
                             color='red', fontsize=8, fontweight='bold')

        # Tooltip annotation (hidden initially)
        self._opt_annot = self.opt_figure.text(
            0.5, 0.01, '', transform=self.opt_figure.transFigure,
            ha='center', va='bottom', fontsize=9, fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.3', fc='#ffffcc', ec='gray', alpha=0.9))
        self._opt_annot.set_visible(False)

        self.opt_figure.tight_layout()
        self.opt_canvas.draw()

        # Connect hover event
        self._opt_motion_cid = self.opt_canvas.mpl_connect(
            'motion_notify_event', self._on_opt_hover)

    def _on_opt_hover(self, event):
        if event.inaxes is None:
            if self._opt_annot.get_visible():
                self._opt_annot.set_visible(False)
                self.opt_canvas.draw_idle()
            return

        x, y = event.xdata, event.ydata

        if event.inaxes is self._opt_ax_left:
            wavelengths, D_scan, mec_map, metric_label, metric_units = self._opt_left_data
            # Find nearest indices
            j = np.searchsorted(wavelengths, x)
            j = np.clip(j, 0, len(wavelengths) - 1)
            i = np.searchsorted(D_scan, y)
            i = np.clip(i, 0, len(D_scan) - 1)
            mec_val = mec_map[i, j]
            if np.isnan(mec_val):
                txt = f'λ={wavelengths[j]:.2f} мкм  D={D_scan[i]:.4f} мкм  {metric_label}=NaN'
            else:
                txt = f'λ={wavelengths[j]:.2f} мкм  D={D_scan[i]:.4f} мкм  {metric_label}={mec_val:.4e} {metric_units}'
        elif event.inaxes is self._opt_ax_right:
            dmax_grid, dmin_grid, window_mec, metric_label, metric_units = self._opt_right_data
            # Log-space search for log-scale axes
            j = np.searchsorted(dmax_grid, x)
            j = np.clip(j, 0, len(dmax_grid) - 1)
            i = np.searchsorted(dmin_grid, y)
            i = np.clip(i, 0, len(dmin_grid) - 1)
            mec_val = window_mec[i, j]
            if np.isnan(mec_val):
                txt = f'D_min={dmin_grid[i]:.4f}  D_max={dmax_grid[j]:.4f}  {metric_label}=NaN'
            else:
                txt = f'D_min={dmin_grid[i]:.4f}  D_max={dmax_grid[j]:.4f}  {metric_label}={mec_val:.4e} {metric_units}'
        else:
            if self._opt_annot.get_visible():
                self._opt_annot.set_visible(False)
                self.opt_canvas.draw_idle()
            return

        self._opt_annot.set_text(txt)
        self._opt_annot.set_visible(True)
        self.opt_canvas.draw_idle()

    def on_optim_finish(self, success, message):
        self.btn_opt_run.setEnabled(True)
        self.btn_opt_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if "Критическая ошибка" in message:
                QMessageBox.critical(self, "Ошибка", message)
            else:
                QMessageBox.warning(self, "Остановка", message)

    def save_optim(self):
        if not self.last_optim_results:
            return
        fn, _ = QFileDialog.getSaveFileName(self, "Экспорт", "optim_results.txt", "Text (*.txt)")
        if not fn:
            return
        try:
            res = self.last_optim_results
            with open(fn, "w", encoding="utf-8") as f:
                f.write("MIE OPTIMIZATION REPORT\n")
                f.write("=======================\n\n")
                f.write(f"D_min* = {res['d_min']:.6f} мкм\n")
                f.write(f"D_max* = {res['d_max']:.6f} мкм\n")
                f.write(f"μ*     = {res['mu']:.6f}\n")
                f.write(f"σ*     = {res['sigma']:.6f}\n")
                f.write(f"L      = {res.get('path_length_m', 1.0):.6f} м\n")
                f.write(f"MEC_mean = {res['mec_mean']:.6e} м²/г\n")
                f.write(f"MEC_min  = {res['mec_min']:.6e} м²/г\n")
                f.write(f"(MEC*L)_mean = {res.get('mec_l_mean', res['mec_mean']):.6e} м³/г\n")
                f.write(f"(MEC*L)_min  = {res.get('mec_l_min', res['mec_min']):.6e} м³/г\n\n")
                f.write(f"{'λ(мкм)':>10} | {'MEC(м²/г)':>15} | {'MEC*L(м³/г)':>15}\n")
                f.write("-" * 50 + "\n")
                mec_l_spectrum = res.get("mec_l_spectrum", res["mec_spectrum"])
                for wl, mec, mec_l in zip(res["wavelengths"], res["mec_spectrum"], mec_l_spectrum):
                    f.write(f"{wl:10.4f} | {mec:15.6e} | {mec_l:15.6e}\n")
            QMessageBox.information(self, "OK", "Сохранено.")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def _spin_ndec(self, val, minv, maxv, dec=1, step=0.1):
        s = QDoubleSpinBox()
        s.setLocale(QLocale(QLocale.C))
        s.setRange(minv, maxv)
        s.setValue(val)
        s.setDecimals(dec)
        s.setSingleStep(step)
        return s

    def _combo_data(self, combo, fallback):
        data = combo.currentData()
        return fallback if data is None else data

    def _on_dist_type_changed(self, _text=None):
        dist_mode = self._combo_data(self.dist_type_combo, DIST_LOGNORMAL)
        is_mono = (dist_mode == DIST_MONODISPERSE)
        is_lognormal = (dist_mode == DIST_LOGNORMAL)
        is_custom = (dist_mode == DIST_CUSTOM)
        is_poly = is_lognormal or is_custom

        self.s_d_mono.setEnabled(is_mono)
        self.lbl_d_mono.setEnabled(is_mono)

        self.s_dmin.setEnabled(is_poly)
        self.s_dmax.setEnabled(is_poly)
        self.s_points.setEnabled(is_poly)
        self.lbl_dmin.setEnabled(is_poly)
        self.lbl_dmax.setEnabled(is_poly)
        self.lbl_points.setEnabled(is_poly)

        self.s_dgm.setEnabled(is_lognormal)
        self.s_sigma.setEnabled(is_lognormal)
        self.lbl_dgm.setEnabled(is_lognormal)
        self.lbl_sigma.setEnabled(is_lognormal)

        self.s_custom_A.setEnabled(is_custom)
        self.s_custom_mu.setEnabled(is_custom)
        self.s_custom_sigma.setEnabled(is_custom)
        self.lbl_custom_A.setEnabled(is_custom)
        self.lbl_custom_mu.setEnabled(is_custom)
        self.lbl_custom_sigma.setEnabled(is_custom)

    def _on_conc_mode_changed(self, _text=None):
        conc_mode = self._combo_data(self.conc_mode, CONC_MASS)
        if conc_mode == CONC_NUMBER:
            self.conc_value.setDecimals(6)
            self.conc_value.setSuffix(" 1/м³")
            if self.conc_value.value() <= 0:
                self.conc_value.setValue(1e9)
        else:
            self.conc_value.setDecimals(2)
            self.conc_value.setSuffix(" г/м³")
            if self.conc_value.value() <= 0:
                self.conc_value.setValue(1.2)

    def _on_inv_wl_mode_changed(self, _text=None):
        wl_mode = self._combo_data(self.inv_wl_mode, INV_WL_SINGLE)
        is_range = (wl_mode == INV_WL_RANGE)
        self.inv_lambda.setEnabled(not is_range)
        self.lbl_inv_lambda.setEnabled(not is_range)
        for w in [
            self.inv_wl_min,
            self.inv_wl_max,
            self.inv_wl_step,
            self.lbl_inv_wl_min,
            self.lbl_inv_wl_max,
            self.lbl_inv_wl_step,
        ]:
            w.setEnabled(is_range)

    def _on_inv_mode_changed(self, _text=None):
        input_mode = self._combo_data(self.inv_input_mode, INV_EFFECTIVE_TRANSMITTANCE)
        needs_mass_conc = inverse_requires_mass_conc(input_mode)
        self.inv_mass_conc.setEnabled(needs_mass_conc)
        self.lbl_inv_mass_conc.setEnabled(needs_mass_conc)
        if inverse_uses_transmittance(input_mode):
            self.inv_target_value.setDecimals(6)
            self.inv_target_value.setSuffix("")
        elif input_mode == INV_TAU:
            self.inv_target_value.setDecimals(6)
            self.inv_target_value.setSuffix("")
        elif input_mode == INV_ALPHA:
            self.inv_target_value.setDecimals(6)
            self.inv_target_value.setSuffix(" 1/м")
        else:
            self.inv_target_value.setDecimals(6)
            self.inv_target_value.setSuffix(" м²/г")

    def log_msg(self, t):
        self.log.append(t)
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    def start(self):
        if self.worker and self.worker.isRunning():
            return

        fracs = {}
        tot = 0.0
        for c, (cb, sp) in self.mats.items():
            if cb.isChecked():
                fracs[c] = float(sp.value())
                tot += float(sp.value())

        if not fracs or tot <= 0:
            QMessageBox.warning(self, "Ошибка", "Не выбран состав смеси.")
            return

        dist_mode = self._combo_data(self.dist_type_combo, DIST_LOGNORMAL)
        is_monodisperse = (dist_mode == DIST_MONODISPERSE)
        is_custom = (dist_mode == DIST_CUSTOM)

        if is_monodisperse:
            d_mono = float(self.s_d_mono.value())
            if d_mono <= 0:
                QMessageBox.warning(self, "Ошибка", "Диаметр должен быть > 0.")
                return
        else:
            dmin, dmax = float(self.s_dmin.value()), float(self.s_dmax.value())
            if dmin >= dmax:
                QMessageBox.warning(self, "Ошибка", "D_min >= D_max.")
                return

        wlmin, wlmax = float(self.s_wl_min.value()), float(self.s_wl_max.value())
        step = float(self.s_wl_step.value())

        if wlmin >= wlmax:
            QMessageBox.warning(self, "Ошибка", "Lambda_min >= Lambda_max.")
            return
        if step <= 0:
            QMessageBox.warning(self, "Ошибка", "Шаг по длине волны должен быть > 0.")
            return

        wls = make_wavelengths(wlmin, wlmax, step)

        if is_monodisperse:
            total_ops = int(len(wls) * len(fracs))
        else:
            N_D = int(self.s_points.value())
            total_ops = int(len(wls) * N_D * len(fracs))

        if total_ops > 2000000:
            res = QMessageBox.question(
                self,
                "Подтверждение",
                f"Расчет потребует ~{total_ops/1e6:.1f} млн вызовов MieQ.\nПродолжить?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                return

        norm_fracs = {k: v / tot for k, v in fracs.items()}

        conc_mode = self._combo_data(self.conc_mode, CONC_MASS)
        conc_value = float(self.conc_value.value())
        if conc_value <= 0 or (not np.isfinite(conc_value)):
            QMessageBox.warning(self, "Ошибка", "Концентрация должна быть > 0.")
            return
        path_length_m = float(self.path_length.value())
        if path_length_m <= 0 or (not np.isfinite(path_length_m)):
            QMessageBox.warning(self, "Ошибка", "Длина трассы L должна быть > 0.")
            return

        if is_monodisperse:
            p = {
                "monodisperse": True,
                "D_um": float(self.s_d_mono.value()),
                "wl_range": (wlmin, wlmax, step),
                "fractions": norm_fracs,
                "conc_mode": conc_mode,
                "conc_value": conc_value,
                "path_length_m": path_length_m,
            }
        elif is_custom:
            p = {
                "monodisperse": False,
                "dist_type": DIST_CUSTOM,
                "d_range": (dmin, dmax),
                "custom_A": float(self.s_custom_A.value()),
                "custom_mu": float(self.s_custom_mu.value()),
                "custom_sigma": float(self.s_custom_sigma.value()),
                "wl_range": (wlmin, wlmax, step),
                "points_d": int(self.s_points.value()),
                "fractions": norm_fracs,
                "conc_mode": conc_mode,
                "conc_value": conc_value,
                "path_length_m": path_length_m,
            }
        else:
            p = {
                "monodisperse": False,
                "dist_type": DIST_LOGNORMAL,
                "d_range": (dmin, dmax),
                "d_dist": (float(self.s_dgm.value()), float(self.s_sigma.value())),
                "wl_range": (wlmin, wlmax, step),
                "points_d": int(self.s_points.value()),
                "fractions": norm_fracs,
                "conc_mode": conc_mode,
                "conc_value": conc_value,
                "path_length_m": path_length_m,
            }

        self.log.clear()
        self.log_msg("=== ЗАПУСК НОВОГО РАСЧЕТА ===")

        self.worker = CalculationWorker(p)
        self.worker.log_signal.connect(self.log_msg)
        self.worker.progress_signal.connect(self.pbar.setValue)
        self.worker.result_signal.connect(self.on_result)
        self.worker.finished_signal.connect(self.on_finish)

        self.btn_run.setEnabled(False)
        self.btn_stop.setEnabled(True)
        self.btn_save.setEnabled(False)
        self.lbl_st.setText("Идет расчет...")
        self.pbar.setValue(0)
        self.worker.start()

    def stop(self):
        if self.worker:
            self.worker.stop()
            self.lbl_st.setText("Прерывание...")

    def on_finish(self, success, message):
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if "Критическая ошибка" in message:
                QMessageBox.critical(self, "Ошибка", message)
            else:
                QMessageBox.warning(self, "Остановка", message)

    def on_result(self, res):
        self.last_results = res
        self.btn_save.setEnabled(True)

    def start_inverse(self):
        if self.inverse_worker and self.inverse_worker.isRunning():
            return

        fracs = {}
        tot = 0.0
        for c, (cb, sp) in self.inv_mats.items():
            if cb.isChecked():
                fracs[c] = float(sp.value())
                tot += float(sp.value())

        if not fracs or tot <= 0:
            QMessageBox.warning(self, "Ошибка", "Не выбран состав смеси.")
            return

        norm_fracs = {k: v / tot for k, v in fracs.items()}

        inv_wl_mode = self._combo_data(self.inv_wl_mode, INV_WL_SINGLE)
        inv_wl_is_range = (inv_wl_mode == INV_WL_RANGE)
        if inv_wl_is_range:
            wl_min = float(self.inv_wl_min.value())
            wl_max = float(self.inv_wl_max.value())
            wl_step = float(self.inv_wl_step.value())
            if wl_min >= wl_max or wl_step <= 0:
                QMessageBox.warning(self, "Ошибка", "Неверные параметры диапазона λ.")
                return
            inv_wavelengths = make_wavelengths(wl_min, wl_max, wl_step)
            if inv_wavelengths.size == 0:
                QMessageBox.warning(self, "Ошибка", "Пустая сетка длин волн.")
                return
            lam_um = None
        else:
            lam_um = float(self.inv_lambda.value())
            if lam_um <= 0:
                QMessageBox.warning(self, "Ошибка", "λ должна быть > 0.")
                return
            wl_min = wl_max = wl_step = None
            inv_wavelengths = np.array([lam_um], dtype=float)

        input_mode = self._combo_data(self.inv_input_mode, INV_EFFECTIVE_TRANSMITTANCE)
        target_value = float(self.inv_target_value.value())
        if target_value <= 0 or (
            inverse_uses_transmittance(input_mode)
            and target_value > 1.0
        ):
            QMessageBox.warning(self, "Ошибка", "Целевое значение должно быть > 0, а для T не больше 1.")
            return
        path_length_m = float(self.inv_path_length.value())
        if path_length_m <= 0 or (not np.isfinite(path_length_m)):
            QMessageBox.warning(self, "Ошибка", "Длина трассы L должна быть > 0.")
            return

        d_min = float(self.inv_d_min.value())
        d_max = float(self.inv_d_max.value())
        if d_min >= d_max:
            QMessageBox.warning(self, "Ошибка", "D_min >= D_max.")
            return

        n_scan = int(self.inv_n_scan.value())
        total_ops = int(n_scan * len(inv_wavelengths))
        if total_ops > 2000000:
            res = QMessageBox.question(
                self,
                "Подтверждение",
                f"Обратная задача потребует ~{total_ops/1e6:.1f} млн расчетов MEC.\nПродолжить?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                return

        p = {
            "fractions": norm_fracs,
            "wl_mode": INV_WL_RANGE if inv_wl_is_range else INV_WL_SINGLE,
            "input_mode": input_mode,
            "target_value": target_value,
            "path_length_m": path_length_m,
            "D_min_um": d_min,
            "D_max_um": d_max,
            "N_scan": n_scan,
        }
        if inv_wl_is_range:
            p["wl_range"] = (wl_min, wl_max, wl_step)
        else:
            p["lambda_um"] = lam_um

        if inverse_requires_mass_conc(input_mode):
            mass_conc_g = float(self.inv_mass_conc.value())
            if mass_conc_g <= 0:
                QMessageBox.warning(self, "Ошибка", "Массовая концентрация должна быть > 0.")
                return
            p["mass_conc_g"] = mass_conc_g

        self.log.clear()

        self.inverse_worker = InverseWorker(p)
        self.inverse_worker.log_signal.connect(self.log_msg)
        self.inverse_worker.progress_signal.connect(self.pbar.setValue)
        self.inverse_worker.result_signal.connect(self.on_inverse_result)
        self.inverse_worker.finished_signal.connect(self.on_inverse_finish)

        self.btn_inv_run.setEnabled(False)
        self.btn_inv_stop.setEnabled(True)
        self.btn_inv_save.setEnabled(False)
        self.lbl_st.setText("Поиск D...")
        self.pbar.setValue(0)
        self.inverse_worker.start()

    def stop_inverse(self):
        if self.inverse_worker:
            self.inverse_worker.stop()
            self.lbl_st.setText("Прерывание...")

    def on_inverse_finish(self, success, message):
        self.btn_inv_run.setEnabled(True)
        self.btn_inv_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if "Критическая ошибка" in message:
                QMessageBox.critical(self, "Ошибка", message)
            else:
                QMessageBox.warning(self, "Остановка", message)

    def on_inverse_result(self, res):
        self.last_inverse_results = res
        self.btn_inv_save.setEnabled(True)

    def save(self):
        if not self.last_results:
            return
        fn, _ = QFileDialog.getSaveFileName(self, "Экспорт", "results.txt", "Text (*.txt)")
        if not fn:
            return
        try:
            p = self.last_results["params"]
            data = self.last_results["data"]
            num_conc = float(self.last_results.get("num_conc", float("nan")))
            mass_conc_g = float(self.last_results.get("mass_conc_g", float("nan")))
            avg_mass = float(self.last_results.get("avg_mass_kg", float("nan")))
            path_length_m = float(p.get("path_length_m", 1.0))

            mats = list(p["fractions"].keys())
            is_mono = p.get("monodisperse", False)

            with open(fn, "w", encoding="utf-8") as f:
                f.write("MIE CALCULATION REPORT (ALPHA + MEC)\n")
                f.write("====================================\n")

                if is_mono:
                    f.write(f"Monodisperse: D={p['D_um']:.4f} um\n")
                elif p.get("dist_type") == DIST_CUSTOM:
                    f.write(f"Dist: Custom(A={p['custom_A']:.6f}, mu={p['custom_mu']:.6f}, sigma={p['custom_sigma']:.6f}), Range=[{p['d_range'][0]}, {p['d_range'][1]}] um, N_D={p['points_d']}\n")
                else:
                    f.write(f"Dist: Log-Normal(Dg={p['d_dist'][0]}, Sig={p['d_dist'][1]}), Range=[{p['d_range'][0]}, {p['d_range'][1]}] um, N_D={p['points_d']}\n")

                f.write(f"Spectrum: [{p['wl_range'][0]}, {p['wl_range'][1]}] um, step={p['wl_range'][2]} um\n")
                f.write(f"Path length L: {path_length_m:.6f} m\n")
                f.write("Mix (Number Fraction):\n")
                for m, v in p["fractions"].items():
                    f.write(f"  {m}: {v*100:.4f}%\n")
                f.write(f"\nAvg particle mass (mixture): {avg_mass:.8e} kg/particle\n")
                f.write(f"Number concentration N: {num_conc:.8e} 1/m^3\n")
                f.write(f"Mass concentration: {mass_conc_g:.8e} g/m^3\n\n")

                if is_mono:
                    hdr = f"{'WL(um)':>10} | {'Cext(um^2)':>15} | {'alpha(1/m)':>15} | {'tau':>15} | {'T':>12} | {'MEC(m^2/g)':>15}"
                else:
                    hdr = f"{'WL(um)':>10} | {'<Cext>(um^2)':>15} | {'alpha(1/m)':>15} | {'tau':>15} | {'T':>12} | {'MEC(m^2/g)':>15}"

                for m in mats:
                    hdr += f" | {('Cext_'+m+'(um^2)'):>15}"
                f.write(hdr + "\n")
                f.write("-" * len(hdr) + "\n")

                for r in data:
                    line = (
                        f"{r['wl']:10.4f} | {r['cext_um2']:15.6e} | {r['alpha_1m']:15.6e} | "
                        f"{r['tau']:15.6e} | {r['transmittance']:12.6e} | {r['mec_m2g']:15.6e}"
                    )
                    for m in mats:
                        line += f" | {r['parts'][m]:15.6e}"
                    f.write(line + "\n")

                summary = summarize_forward_rows(data, mats)
                if summary:
                    f.write("-" * len(hdr) + "\n")
                    avg_line = (
                        f"{'AVG':>10} | {summary['cext_um2']:15.6e} | {summary['alpha_1m']:15.6e} | "
                        f"{summary['tau']:15.6e} | {summary['transmittance']:12.6e} | {summary['mec_m2g']:15.6e}"
                    )
                    for m in mats:
                        avg_line += f" | {summary['parts'][m]:15.6e}"
                    f.write(avg_line + "\n")
                    f.write(f"T_eff=exp(-AVG tau): {summary['eff_transmittance']:.6e}\n")

            QMessageBox.information(self, "OK", "Сохранено.")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def save_inverse(self):
        if not self.last_inverse_results:
            return
        fn, _ = QFileDialog.getSaveFileName(self, "Экспорт", "inverse_results.txt", "Text (*.txt)")
        if not fn:
            return
        try:
            res = self.last_inverse_results
            p = res["params"]
            solutions = res["solutions"]
            scan_D = res["scan_D"]
            scan_MEC = res["scan_MEC"]
            target_mec = res["target_mec"]
            metric_label = res.get("metric_label", "MEC")
            metric_units = res.get("metric_units", "м²/г")
            wavelengths = res.get("wavelengths", [])
            input_mode = p["input_mode"]
            use_avg_spectrum = p.get("wl_mode") == INV_WL_RANGE
            alpha_name, tau_name = inverse_solution_names(use_avg_spectrum)
            transmittance_note = (
                "reference"
                if inverse_transmittance_is_reference(input_mode, use_avg_spectrum)
                else "check"
            )
            units_suffix = f" {metric_units}" if metric_units else ""

            with open(fn, "w", encoding="utf-8") as f:
                f.write("MIE INVERSE PROBLEM REPORT (MONODISPERSE)\n")
                f.write("=========================================\n\n")

                if use_avg_spectrum:
                    f.write(f"Spectrum: [{p['wl_range'][0]}, {p['wl_range'][1]}] мкм, step={p['wl_range'][2]} мкм, N={len(wavelengths)}\n")
                else:
                    f.write(f"λ = {p['lambda_um']:.4f} мкм\n")
                f.write(f"L = {p.get('path_length_m', 1.0):.6f} м\n")
                f.write(f"Input mode: {input_mode}\n")
                f.write(f"Target value: {p['target_value']:.6e}\n")
                if inverse_requires_mass_conc(input_mode):
                    f.write(f"Mass concentration: {p['mass_conc_g']:.6e} г/м³\n")
                if inverse_uses_transmittance(input_mode) and res.get('equivalent_mec') is not None:
                    f.write(f"Equivalent -ln(T)/(rho_mass*L): {res['equivalent_mec']:.6e} м²/г\n")
                f.write(f"Target {metric_label}: {target_mec:.6e}{units_suffix}\n")
                f.write(f"Search range: [{p['D_min_um']:.4f}, {p['D_max_um']:.4f}] мкм\n")
                f.write(f"Scan points: {p['N_scan']}\n\n")

                f.write("Mix (Number Fraction):\n")
                for m, v in p["fractions"].items():
                    f.write(f"  {m}: {v*100:.4f}%\n")
                f.write("\n")

                f.write(f"SOLUTIONS FOUND: {len(solutions)}\n")
                f.write("-" * 50 + "\n")

                for i, sol in enumerate(solutions):
                    f.write(f"\nSolution #{i+1}:\n")
                    f.write(f"  D = {sol['D_um']:.6f} мкм\n")
                    if sol['N'] is not None:
                        f.write(f"  N = {sol['N']:.6e} 1/м³\n")
                        f.write(f"  {alpha_name} (check) = {sol.get('alpha_check', float('nan')):.6e} 1/м\n")
                        f.write(f"  {tau_name}=alpha*L (check) = {sol.get('tau_check', float('nan')):.6e}\n")
                        f.write(f"  exp(-{tau_name}) ({transmittance_note}) = {sol.get('transmittance_check', float('nan')):.6e}\n")
                    if inverse_shows_avg_mec_reference(input_mode):
                        f.write(f"  AVG MEC (reference) = {sol.get('avg_mec_check', float('nan')):.6e} м²/г\n")
                    f.write(f"  {metric_label} (check) = {sol['mec_check']:.6e}{units_suffix}\n")
                    f.write(f"  Relative error = {sol.get('rel_error', 0)*100:.4f}%\n")
                    f.write(f"  m_particle = {sol['m_particle_kg']:.6e} кг\n")

                f.write(f"\n\nSCAN DATA (D vs {metric_label}):\n")
                f.write("-" * 50 + "\n")
                metric_header = metric_label if not metric_units else f"{metric_label}({metric_units})"
                f.write(f"{'D(um)':>12} | {metric_header:>18}\n")
                f.write("-" * 30 + "\n")

                step = max(1, len(scan_D) // 100)
                for i in range(0, len(scan_D), step):
                    f.write(f"{scan_D[i]:12.6f} | {scan_MEC[i]:18.6e}\n")

            QMessageBox.information(self, "OK", "Сохранено.")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    w = MainWindow()
    w.show()
    sys.exit(app.exec())
