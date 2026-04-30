"""Pure-Python optical / Mie computational core for mie_shield.

No Qt or matplotlib dependencies; this module is safe to import in tests and
in any non-GUI context. The Qt-based workers, GUI, and main entry point live
in mie_shield.py.

Public surface (consumed by mie_shield.py and tests):

* MATERIALS_DB
* TypeAliases: ConcMode, DistributionMode, InverseInputMode,
  InverseWavelengthMode, OptimizationMode, OptimizationCriterion,
  RIModel, Fractions, DensityMap, ForwardRow
* Mode constants: CONC_MASS/CONC_NUMBER, DIST_*, INV_WL_*, INV_*,
  OPT_*, plus DENSITY_FALLBACK and TAU_UNDERFLOW_LIMIT
* Core functions: make_wavelengths, get_ri, material_density_map,
  mixture_density, monodisperse_particle_mass_kg,
  distributed_particle_mass_kg, resolve_concentration,
  transmittance_from_tau, make_forward_row, summarize_forward_rows,
  inverse_requires_mass_conc, inverse_uses_transmittance,
  inverse_metric_label, inverse_metric_units, inverse_solution_names,
  inverse_transmittance_is_reference, inverse_shows_avg_mec_reference,
  resolve_inverse_target, inverse_metric_from_mec_values,
  lognormal_pdf, custom_pdf, safe_mie_qext, qext_to_cext_um2,
  compute_qext_avg, compute_mec_for_d, nan_safe_bisect
"""

from typing import Callable, Iterable, Literal, TypeAlias

import numpy as np
import scipy.integrate

if not hasattr(scipy.integrate, "trapz"):
    scipy.integrate.trapz = scipy.integrate.trapezoid

from PyMieScatt import MieQ

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

ConcMode: TypeAlias = Literal["mass", "number"]
DistributionMode: TypeAlias = Literal["monodisperse", "lognormal", "custom"]
InverseWavelengthMode: TypeAlias = Literal["single", "range"]
InverseInputMode: TypeAlias = Literal["mec", "alpha", "tau", "transmittance", "effective_transmittance"]
OptimizationMode: TypeAlias = Literal["window_only", "full"]
OptimizationCriterion: TypeAlias = Literal["mean", "min"]
RIModel: TypeAlias = Callable[[float], complex]
Fractions: TypeAlias = dict[str, float]
DensityMap: TypeAlias = dict[str, float]
ForwardRow: TypeAlias = dict[str, object]

CONC_MASS: ConcMode = "mass"
CONC_NUMBER: ConcMode = "number"

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
    if conc_mode not in (CONC_MASS, CONC_NUMBER):
        raise ValueError(f"Неизвестный режим концентрации: {conc_mode!r}.")
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
