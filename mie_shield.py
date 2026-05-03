import sys
import traceback
import types
from pathlib import Path

import numpy as np


def _install_scipy_stats_placeholder() -> None:
    # Optional SciPy optimize/interpolate branches import scipy.stats; this app
    # does not use those branches, and excluding stats keeps Nuitka builds sane.
    if "scipy.stats" not in sys.modules:
        sys.modules["scipy.stats"] = types.ModuleType("scipy.stats")


_install_scipy_stats_placeholder()

import scipy.optimize

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PySide6.QtCore import Qt, QThread, QLocale, Signal
from PySide6.QtGui import QAction, QActionGroup, QFont, QIcon
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

APP_ID = "com.local.MieShield"
APP_NAME = "MieShield"

from mie_i18n import i18n, t
from mie_i18n_strings import LANGUAGE_NAMES, SUPPORTED_LANGUAGES


def material_label(code: str) -> str:
    return f"{code} - {t(f'material.{code}')}"


def _resource_path(name: str) -> Path:
    return Path(__file__).resolve().with_name(name)


def _application_icon() -> QIcon:
    icon_path = _resource_path("icon.png")
    if icon_path.exists():
        return QIcon(str(icon_path))
    return QIcon()


def _set_windows_app_user_model_id() -> None:
    if sys.platform != "win32":
        return
    try:
        import ctypes

        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(APP_ID)
    except Exception:
        pass

# Explicit core imports keep the boundary between numerical code and the
# Qt worker/UI layer visible.
from mie_core import (
    # Reference data
    MATERIALS_DB,
    MieCoreError,
    # Mode constants (typed string literals)
    CONC_MASS,
    CONC_NUMBER,
    DIST_MONODISPERSE,
    DIST_LOGNORMAL,
    DIST_CUSTOM,
    INV_WL_SINGLE,
    INV_WL_RANGE,
    INV_MEC,
    INV_ALPHA,
    INV_TAU,
    INV_TRANSMITTANCE,
    INV_EFFECTIVE_TRANSMITTANCE,
    OPT_WINDOW_ONLY,
    OPT_FULL,
    OPT_MEAN,
    OPT_MIN,
    # Spectrum / RI
    make_wavelengths,
    get_ri,
    # Density / mass helpers
    material_density_map,
    mixture_density,
    monodisperse_particle_mass_kg,
    distributed_particle_mass_kg,
    # Concentration / forward row helpers
    resolve_concentration,
    transmittance_from_tau,
    make_forward_row,
    summarize_forward_rows,
    # Inverse-problem helpers
    inverse_requires_mass_conc,
    inverse_uses_transmittance,
    inverse_metric_label,
    inverse_metric_units,
    inverse_solution_names,
    inverse_transmittance_is_reference,
    inverse_shows_avg_mec_reference,
    resolve_inverse_target,
    inverse_metric_from_mec_values,
    # Distributions
    lognormal_pdf,
    lognormal_cdf,
    custom_pdf,
    # Mie call wrappers
    safe_mie_qext,
    qext_to_cext_um2,
    compute_qext_avg,
    compute_mec_for_d,
    # Numerical helpers
    nan_safe_bisect,
    trapz,
)

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


class CalculationWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str, object)

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
                raise MieCoreError("err.empty_wavelength_grid")

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
                raise MieCoreError("err.path_length_positive")
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
                        self.finished_signal.emit(False, "calc.status_aborted", {})
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
                        self.finished_signal.emit(
                            False,
                            "calc.err_mie_mono",
                            {"lam_um": lam_um, "lost_pct": mixture_lost_weight * 100},
                        )
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
                self.finished_signal.emit(True, "calc.status_success", {})

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
                    mass_numerical = trapz(pdf_values, diameters_um)
                    if not np.isfinite(mass_numerical) or mass_numerical <= 1e-12:
                        raise MieCoreError("err.pdf_integral_invalid")

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
                    mass_numerical = trapz(pdf_values, diameters_um)
                    if not np.isfinite(mass_numerical) or mass_numerical <= 1e-12:
                        raise MieCoreError("err.pdf_integral_invalid")

                    shape_param = np.log(sigma_g)
                    cdf_max = lognormal_cdf(d_max_um, dg_um, sigma_g)
                    cdf_min = lognormal_cdf(d_min_um, dg_um, sigma_g)
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
                        self.finished_signal.emit(False, "calc.status_aborted", {})
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
                        avg_c_ext_mat = float(trapz(integrand, diameters_um))

                        c_ext_parts[mat_code] = avg_c_ext_mat
                        c_ext_total_mix += float(num_fraction) * avg_c_ext_mat

                        if np.any(mat_failed_mask):
                            pdf_failed = np.where(mat_failed_mask, pdf_normalized, 0.0)
                            w_lost_mat = float(trapz(pdf_failed, diameters_um))
                            mixture_lost_weight += float(num_fraction) * w_lost_mat

                    if mixture_lost_weight > 0.02:
                        self.finished_signal.emit(
                            False,
                            "calc.err_mie_distribution",
                            {"lam_um": lam_um, "lost_pct": mixture_lost_weight * 100},
                        )
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
                self.finished_signal.emit(True, "calc.status_success", {})

        except MieCoreError as exc:
            self.finished_signal.emit(False, exc.code, exc.params)
        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, "err.critical", {"traceback": err_msg})


class InverseWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str, object)

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
                    raise MieCoreError("err.empty_wavelength_grid")
                actual_wl_step = (wl_max - wl_min) / (wavelengths.size - 1) if wavelengths.size > 1 else wl_step
            else:
                lam_um = float(p["lambda_um"])
                wavelengths = np.array([lam_um], dtype=float)
                actual_wl_step = 0.0
            input_mode = p["input_mode"]
            target_value = float(p["target_value"])
            path_length_m = float(p.get("path_length_m", 1.0))
            if not np.isfinite(path_length_m) or path_length_m <= 0:
                raise MieCoreError("err.path_length_positive")

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
                    self.finished_signal.emit(False, "calc.status_aborted", {})
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
                    self.finished_signal.emit(False, "calc.status_aborted", {})
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
            self.finished_signal.emit(True, "inverse.status_finished", {"count": len(solutions)})

        except MieCoreError as exc:
            self.finished_signal.emit(False, exc.code, exc.params)
        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, "err.critical", {"traceback": err_msg})


class OptimizationWorker(QThread):
    log_signal = Signal(str)
    progress_signal = Signal(int)
    result_signal = Signal(dict)
    finished_signal = Signal(bool, str, object)

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
                self.finished_signal.emit(False, "err.empty_wavelength_grid", {})
                return
            path_length_m = float(p.get("path_length_m", 1.0))
            if not np.isfinite(path_length_m) or path_length_m <= 0:
                raise MieCoreError("err.path_length_positive")

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
                    self.finished_signal.emit(False, "calc.status_aborted", {})
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
                    self.finished_signal.emit(False, "calc.status_aborted", {})
                    return
                for j, dx in enumerate(dmax_grid):
                    done_w += 1
                    if dx <= dm * 1.01:
                        continue
                    d_g = np.geomspace(dm, dx, N_D_points)
                    pv = custom_pdf(d_g, 1.0, mu_for_map, sigma_for_map)
                    ps = trapz(pv, d_g)
                    if ps <= 1e-30:
                        continue
                    pn = pv / ps
                    log_d_g = np.log(d_g)
                    mec_wl = np.zeros(n_wl)
                    for k in range(n_wl):
                        pts = np.column_stack([log_d_g, np.full(N_D_points, wavelengths[k])])
                        mec_wl[k] = trapz(interp(pts) * pn, d_g)
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
                pdf_sum = trapz(pdf_vals, d_grid)
                if pdf_sum <= 1e-30:
                    return 1e10
                pdf_norm = pdf_vals / pdf_sum

                log_d_grid = np.log(d_grid)
                mec_per_wl = np.zeros(n_wl)
                for j in range(n_wl):
                    pts = np.column_stack([log_d_grid, np.full(N_D_points, wavelengths[j])])
                    mec_mono = interp(pts)
                    mec_per_wl[j] = trapz(mec_mono * pdf_norm, d_grid)

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
                init="random",
                polish=True, disp=False,
            )

            if self.is_aborted:
                self.finished_signal.emit(False, "calc.status_aborted", {})
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
            pdf_sum_best = trapz(pdf_best, d_grid_best)
            pdf_norm_best = pdf_best / pdf_sum_best if pdf_sum_best > 1e-30 else pdf_best

            mec_spectrum = np.zeros(n_wl)
            log_d_best = np.log(d_grid_best)
            for j in range(n_wl):
                pts = np.column_stack([log_d_best, np.full(N_D_points, wavelengths[j])])
                mec_mono = interp(pts)
                mec_spectrum[j] = trapz(mec_mono * pdf_norm_best, d_grid_best)

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
            self.finished_signal.emit(
                True,
                "optim.status_finished",
                {"mec_l_mean": mec_l_mean, "mec_l_min": mec_l_min},
            )

        except MieCoreError as exc:
            self.finished_signal.emit(False, exc.code, exc.params)
        except Exception:
            err_msg = traceback.format_exc()
            self.finished_signal.emit(False, "err.critical", {"traceback": err_msg})


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(t("app.title"))
        self.setWindowIcon(_application_icon())
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
        self._language_menu = None
        self._language_actions = {}
        self._material_checkboxes = []
        self._text_bindings = []
        self._combo_bindings = []

        self._build_language_menu()

        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        left = QVBoxLayout()
        root.addLayout(left, 1)
        right = QVBoxLayout()
        root.addLayout(right, 2)

        self.tabs = QTabWidget()
        left.addWidget(self.tabs)

        self.tab_forward = QWidget()
        self.tabs.addTab(self.tab_forward, "")
        self._build_forward_tab(self.tab_forward)

        self.tab_inverse = QWidget()
        self.tabs.addTab(self.tab_inverse, "")
        self._build_inverse_tab(self.tab_inverse)

        self.tab_optim = QWidget()
        self.tabs.addTab(self.tab_optim, "")
        self._build_optim_tab(self.tab_optim)

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

        self.lbl_st = QLabel("")
        right.addWidget(self.lbl_st)

        self.tabs.currentChanged.connect(self._on_tab_changed)
        i18n.language_changed.connect(self._retranslate)
        self._retranslate()

    def _build_language_menu(self):
        self._language_menu = self.menuBar().addMenu("")
        group = QActionGroup(self)
        group.setExclusive(True)
        for code in SUPPORTED_LANGUAGES:
            action = QAction(LANGUAGE_NAMES[code], self)
            action.setCheckable(True)
            action.triggered.connect(lambda _checked=False, lang=code: i18n.set_language(lang))
            group.addAction(action)
            self._language_menu.addAction(action)
            self._language_actions[code] = action

    def _retranslate(self):
        self.setWindowTitle(t("app.title"))
        if self._language_menu is not None:
            self._language_menu.setTitle(t("menu.language"))
        for code, action in self._language_actions.items():
            action.setText(LANGUAGE_NAMES[code])
            action.setChecked(code == i18n.language)
        for widget, method_name, key in self._text_bindings:
            getattr(widget, method_name)(t(key))
        for combo, items in self._combo_bindings:
            current_data = combo.currentData()
            was_blocked = combo.blockSignals(True)
            for index, (key, _data) in enumerate(items):
                combo.setItemText(index, t(key))
            index = combo.findData(current_data)
            if index >= 0:
                combo.setCurrentIndex(index)
            combo.blockSignals(was_blocked)
        self.tabs.setTabText(self.tabs.indexOf(self.tab_forward), t("tab.forward"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_inverse), t("tab.inverse"))
        self.tabs.setTabText(self.tabs.indexOf(self.tab_optim), t("tab.optim"))
        for code, checkbox in self._material_checkboxes:
            checkbox.setText(material_label(code))
        if hasattr(self, "conc_mode"):
            self._on_conc_mode_changed()
        if hasattr(self, "inv_input_mode"):
            self._on_inv_mode_changed()
        self.lbl_st.setText(t("status.ready"))

    def _message_from_key(self, key: str, params: object | None = None) -> str:
        kwargs = params if isinstance(params, dict) else {}
        return t(key, **kwargs)

    def _bind_text(self, widget, method_name: str, key: str):
        self._text_bindings.append((widget, method_name, key))
        getattr(widget, method_name)(t(key))
        return widget

    def _tr_label(self, key: str) -> QLabel:
        return self._bind_text(QLabel(""), "setText", key)

    def _tr_groupbox(self, key: str) -> QGroupBox:
        return self._bind_text(QGroupBox(""), "setTitle", key)

    def _tr_button(self, key: str) -> QPushButton:
        return self._bind_text(QPushButton(""), "setText", key)

    def _tr_suffix(self, widget, key: str):
        return self._bind_text(widget, "setSuffix", key)

    def _add_tr_row(self, layout: QFormLayout, key: str, widget) -> None:
        layout.addRow(self._tr_label(key), widget)

    def _add_combo_items(self, combo: QComboBox, items: list[tuple[str, object]]) -> None:
        for _key, data in items:
            combo.addItem("", data)
        self._combo_bindings.append((combo, items))

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

        g1 = self._tr_groupbox("forward.group.particles")
        layout.addWidget(g1)
        l1 = QFormLayout(g1)

        self.dist_type_combo = QComboBox()
        self._add_combo_items(
            self.dist_type_combo,
            [
                ("dist.monodisperse", DIST_MONODISPERSE),
                ("dist.lognormal", DIST_LOGNORMAL),
                ("dist.custom", DIST_CUSTOM),
            ],
        )
        self.dist_type_combo.setCurrentIndex(self.dist_type_combo.findData(DIST_LOGNORMAL))
        self.dist_type_combo.currentTextChanged.connect(self._on_dist_type_changed)
        self._add_tr_row(l1, "label.distribution_type", self.dist_type_combo)

        self.s_d_mono = self._spin_ndec(DEFAULTS["D_MONO"], 0.001, 100, dec=3, step=0.01)
        self.s_d_mono.setEnabled(False)
        self.lbl_d_mono = self._tr_label("label.d_um")
        self.lbl_d_mono.setEnabled(False)
        l1.addRow(self.lbl_d_mono, self.s_d_mono)

        self.s_dmin = self._spin_ndec(DEFAULTS["D_MIN"], 0.001, 100, dec=2, step=0.01)
        self.s_dmax = self._spin_ndec(DEFAULTS["D_MAX"], 0.001, 100, dec=2, step=0.01)
        self.s_dgm = self._spin_ndec(DEFAULTS["D_GEOMETRIC_MEAN"], 0.001, 100, dec=2, step=0.01)
        self.s_sigma = self._spin_ndec(DEFAULTS["SIGMA_GEOMETRIC"], 1.01, 10, dec=2, step=0.01)
        self.s_points = QSpinBox()
        self.s_points.setRange(50, 10000)
        self.s_points.setValue(DEFAULTS["POINTS_D"])

        self.lbl_dmin = self._tr_label("label.d_min_um")
        self.lbl_dmax = self._tr_label("label.d_max_um")
        self.lbl_dgm = self._tr_label("label.dg_um")
        self.lbl_sigma = self._tr_label("label.sigma_g")
        self.lbl_points = self._tr_label("label.grid_points")

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

        self.lbl_custom_A = self._tr_label("label.custom_a")
        self.lbl_custom_mu = self._tr_label("label.custom_mu")
        self.lbl_custom_sigma = self._tr_label("label.custom_sigma")
        self.lbl_custom_A.setEnabled(False)
        self.lbl_custom_mu.setEnabled(False)
        self.lbl_custom_sigma.setEnabled(False)

        l1.addRow(self.lbl_custom_A, self.s_custom_A)
        l1.addRow(self.lbl_custom_mu, self.s_custom_mu)
        l1.addRow(self.lbl_custom_sigma, self.s_custom_sigma)

        l1.addRow(self.lbl_points, self.s_points)

        g_wl = self._tr_groupbox("forward.group.spectrum")
        layout.addWidget(g_wl)
        l_wl = QFormLayout(g_wl)
        self.s_wl_min = self._spin_ndec(DEFAULTS["WAVELENGTH_MIN"], 0.1, 1000, dec=2, step=0.01)
        self.s_wl_max = self._spin_ndec(DEFAULTS["WAVELENGTH_MAX"], 0.1, 1000, dec=2, step=0.01)
        self.s_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.01, 100, dec=2, step=0.01)
        self._add_tr_row(l_wl, "label.lambda_min_um", self.s_wl_min)
        self._add_tr_row(l_wl, "label.lambda_max_um", self.s_wl_max)
        self._add_tr_row(l_wl, "label.step_um", self.s_wl_step)

        g2 = self._tr_groupbox("group.mixture_number_fraction")
        layout.addWidget(g2)
        l2 = QVBoxLayout(g2)
        self.mats = {}
        for c, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox("")
            self._material_checkboxes.append((c, cb))
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

        g_conc = self._tr_groupbox("forward.group.concentration")
        layout.addWidget(g_conc)
        l_conc = QFormLayout(g_conc)
        self.conc_mode = QComboBox()
        self._add_combo_items(
            self.conc_mode,
            [
                ("conc.number", CONC_NUMBER),
                ("conc.mass", CONC_MASS),
            ],
        )
        self.conc_mode.setCurrentIndex(self.conc_mode.findData(CONC_MASS))
        self.conc_value = QDoubleSpinBox()
        self.conc_value.setLocale(QLocale(QLocale.C))
        self.conc_value.setRange(1e-30, 1e18)
        self.conc_value.setDecimals(2)
        self.conc_value.setValue(1.2)
        self.conc_value.setAlignment(Qt.AlignRight)
        self._tr_suffix(self.conc_value, "unit.g_per_m3")
        self.path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.1, 1000000, dec=1, step=0.1)
        self._tr_suffix(self.path_length, "unit.m")
        self._add_tr_row(l_conc, "label.type", self.conc_mode)
        self._add_tr_row(l_conc, "label.value", self.conc_value)
        self._add_tr_row(l_conc, "label.path_length", self.path_length)
        self.conc_mode.currentTextChanged.connect(self._on_conc_mode_changed)

        g3 = self._tr_groupbox("forward.group.controls")
        layout.addWidget(g3)
        l3 = QVBoxLayout(g3)
        self.btn_run = self._tr_button("button.run")
        self.btn_run.setFixedHeight(40)
        self.btn_run.clicked.connect(self.start)
        self.btn_stop = self._tr_button("button.stop")
        self.btn_stop.setEnabled(False)
        self.btn_stop.clicked.connect(self.stop)
        self.btn_save = self._tr_button("button.export_result")
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save)
        l3.addWidget(self.btn_run)
        l3.addWidget(self.btn_stop)
        l3.addWidget(self.btn_save)

        layout.addStretch()

    def _build_inverse_tab(self, tab):
        layout = QVBoxLayout(tab)

        g_input = self._tr_groupbox("inverse.group.input")
        layout.addWidget(g_input)
        l_input = QFormLayout(g_input)

        self.inv_wl_mode = QComboBox()
        self._add_combo_items(
            self.inv_wl_mode,
            [
                ("wl.single", INV_WL_SINGLE),
                ("wl.range", INV_WL_RANGE),
            ],
        )
        self.inv_wl_mode.currentTextChanged.connect(self._on_inv_wl_mode_changed)
        self._add_tr_row(l_input, "label.spectrum", self.inv_wl_mode)

        self.inv_lambda = self._spin_ndec(0.59, 0.1, 1000, dec=2, step=0.01)
        self.lbl_inv_lambda = self._tr_label("label.lambda_um")
        l_input.addRow(self.lbl_inv_lambda, self.inv_lambda)

        self.inv_wl_min = self._spin_ndec(0.4, 0.1, 1000, dec=2, step=0.01)
        self.inv_wl_max = self._spin_ndec(0.78, 0.1, 1000, dec=2, step=0.01)
        self.inv_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.01, 100, dec=2, step=0.01)
        self.lbl_inv_wl_min = self._tr_label("label.lambda_min_um_short")
        self.lbl_inv_wl_max = self._tr_label("label.lambda_max_um_short")
        self.lbl_inv_wl_step = self._tr_label("label.lambda_step_um")
        l_input.addRow(self.lbl_inv_wl_min, self.inv_wl_min)
        l_input.addRow(self.lbl_inv_wl_max, self.inv_wl_max)
        l_input.addRow(self.lbl_inv_wl_step, self.inv_wl_step)
        self._on_inv_wl_mode_changed(self.inv_wl_mode.currentText())

        self.inv_path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.1, 1000000, dec=1, step=0.1)
        self._tr_suffix(self.inv_path_length, "unit.m")
        self._add_tr_row(l_input, "label.path_length", self.inv_path_length)

        self.inv_input_mode = QComboBox()
        self._add_combo_items(
            self.inv_input_mode,
            [
                ("inverse.input.effective_transmittance", INV_EFFECTIVE_TRANSMITTANCE),
                ("inverse.input.transmittance", INV_TRANSMITTANCE),
                ("inverse.input.tau", INV_TAU),
                ("inverse.input.alpha", INV_ALPHA),
                ("inverse.input.mec", INV_MEC),
            ],
        )
        self.inv_input_mode.currentTextChanged.connect(self._on_inv_mode_changed)
        self._add_tr_row(l_input, "label.input_type", self.inv_input_mode)

        self.inv_target_value = QDoubleSpinBox()
        self.inv_target_value.setLocale(QLocale(QLocale.C))
        self.inv_target_value.setRange(1e-30, 1e30)
        self.inv_target_value.setDecimals(6)
        self.inv_target_value.setValue(1e-3)
        self.inv_target_value.setAlignment(Qt.AlignRight)
        self._add_tr_row(l_input, "label.value", self.inv_target_value)

        self.inv_mass_conc = self._spin_ndec(1.2, 1e-30, 1e18, dec=6, step=0.001)
        self._tr_suffix(self.inv_mass_conc, "unit.g_per_m3")
        self.lbl_inv_mass_conc = self._tr_label("label.mass_conc")
        l_input.addRow(self.lbl_inv_mass_conc, self.inv_mass_conc)
        self._on_inv_mode_changed(self.inv_input_mode.currentText())

        g_search = self._tr_groupbox("inverse.group.search_range")
        layout.addWidget(g_search)
        l_search = QFormLayout(g_search)

        self.inv_d_min = self._spin_ndec(DEFAULTS["INV_D_MIN"], 0.001, 1000, dec=3, step=0.01)
        self.inv_d_max = self._spin_ndec(DEFAULTS["INV_D_MAX"], 0.001, 1000, dec=3, step=1.0)
        self.inv_n_scan = QSpinBox()
        self.inv_n_scan.setRange(100, 100000)
        self.inv_n_scan.setValue(DEFAULTS["INV_N_SCAN"])

        self._add_tr_row(l_search, "label.d_min_um", self.inv_d_min)
        self._add_tr_row(l_search, "label.d_max_um", self.inv_d_max)
        self._add_tr_row(l_search, "label.scan_points", self.inv_n_scan)

        g_mat = self._tr_groupbox("group.mixture_number_fraction")
        layout.addWidget(g_mat)
        l_mat = QVBoxLayout(g_mat)
        self.inv_mats = {}
        for c, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox("")
            self._material_checkboxes.append((c, cb))
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

        g_ctrl = self._tr_groupbox("inverse.group.controls")
        layout.addWidget(g_ctrl)
        l_ctrl = QVBoxLayout(g_ctrl)

        self.btn_inv_run = self._tr_button("button.search_d")
        self.btn_inv_run.setFixedHeight(40)
        self.btn_inv_run.clicked.connect(self.start_inverse)
        self.btn_inv_stop = self._tr_button("button.stop")
        self.btn_inv_stop.setEnabled(False)
        self.btn_inv_stop.clicked.connect(self.stop_inverse)
        self.btn_inv_save = self._tr_button("button.export_result")
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

        g_wl = self._tr_groupbox("optim.group.spectrum_path_materials")
        layout.addWidget(g_wl)
        l_wl = QFormLayout(g_wl)
        self.opt_wl_min = self._spin_ndec(DEFAULTS["WAVELENGTH_MIN"], 0.1, 1000, dec=2, step=0.01)
        self.opt_wl_max = self._spin_ndec(DEFAULTS["WAVELENGTH_MAX"], 0.1, 1000, dec=2, step=0.01)
        self.opt_wl_step = self._spin_ndec(DEFAULTS["WAVELENGTH_STEP"], 0.01, 100, dec=2, step=0.01)
        self._add_tr_row(l_wl, "label.lambda_min_um_short", self.opt_wl_min)
        self._add_tr_row(l_wl, "label.lambda_max_um_short", self.opt_wl_max)
        self._add_tr_row(l_wl, "label.lambda_step_um", self.opt_wl_step)
        self.opt_path_length = self._spin_ndec(DEFAULTS["PATH_LENGTH_M"], 0.1, 1000000, dec=1, step=0.1)
        self._tr_suffix(self.opt_path_length, "unit.m")
        self._add_tr_row(l_wl, "label.path_length", self.opt_path_length)

        g_mat = self._tr_groupbox("optim.group.mixture_number_fraction")
        layout.addWidget(g_mat)
        l_mat = QVBoxLayout(g_mat)
        self.opt_mats = {}
        for c, _rho in MATERIALS_DB:
            h = QHBoxLayout()
            cb = QCheckBox("")
            self._material_checkboxes.append((c, cb))
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

        g_dist = self._tr_groupbox("optim.group.distribution_params")
        layout.addWidget(g_dist)
        l_dist = QFormLayout(g_dist)
        self.opt_mode = QComboBox()
        self._add_combo_items(
            self.opt_mode,
            [
                ("optim.mode.window_only", OPT_WINDOW_ONLY),
                ("optim.mode.full", OPT_FULL),
            ],
        )
        self.opt_mode.currentTextChanged.connect(self._on_optim_mode_changed)
        self._add_tr_row(l_dist, "label.mode", self.opt_mode)

        self.opt_mu = self._spin_ndec(DEFAULTS["CUSTOM_MU"], -10, 10, dec=6, step=0.01)
        self.opt_sigma = self._spin_ndec(DEFAULTS["CUSTOM_SIGMA"], 0.001, 10, dec=6, step=0.01)
        self.lbl_opt_mu = self._tr_label("label.mu_fixed")
        self.lbl_opt_sigma = self._tr_label("label.sigma_fixed")
        l_dist.addRow(self.lbl_opt_mu, self.opt_mu)
        l_dist.addRow(self.lbl_opt_sigma, self.opt_sigma)

        self.opt_mu_min = self._spin_ndec(-2.0, -10, 10, dec=6, step=0.1)
        self.opt_mu_max = self._spin_ndec(3.0, -10, 10, dec=6, step=0.1)
        self.opt_sigma_min = self._spin_ndec(0.1, 0.001, 10, dec=6, step=0.01)
        self.opt_sigma_max = self._spin_ndec(2.0, 0.001, 10, dec=6, step=0.01)
        self.lbl_opt_mu_min = self._tr_label("label.mu_min")
        self.lbl_opt_mu_max = self._tr_label("label.mu_max")
        self.lbl_opt_sigma_min = self._tr_label("label.sigma_min")
        self.lbl_opt_sigma_max = self._tr_label("label.sigma_max")
        l_dist.addRow(self.lbl_opt_mu_min, self.opt_mu_min)
        l_dist.addRow(self.lbl_opt_mu_max, self.opt_mu_max)
        l_dist.addRow(self.lbl_opt_sigma_min, self.opt_sigma_min)
        l_dist.addRow(self.lbl_opt_sigma_max, self.opt_sigma_max)
        for w in [self.opt_mu_min, self.opt_mu_max, self.opt_sigma_min, self.opt_sigma_max,
                  self.lbl_opt_mu_min, self.lbl_opt_mu_max, self.lbl_opt_sigma_min, self.lbl_opt_sigma_max]:
            w.setEnabled(False)

        g_search = self._tr_groupbox("optim.group.d_search")
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
        self._add_tr_row(l_search, "label.d_scan_min_um", self.opt_d_scan_min)
        self._add_tr_row(l_search, "label.d_scan_max_um", self.opt_d_scan_max)
        self._add_tr_row(l_search, "label.map_points", self.opt_n_d_scan)
        self._add_tr_row(l_search, "label.pdf_points", self.opt_n_d_points)
        self._add_tr_row(l_search, "label.window_grid", self.opt_n_window)

        g_crit = self._tr_groupbox("optim.group.criterion")
        layout.addWidget(g_crit)
        l_crit = QFormLayout(g_crit)
        self.opt_criterion = QComboBox()
        self._add_combo_items(
            self.opt_criterion,
            [
                ("optim.criterion.mean", OPT_MEAN),
                ("optim.criterion.min", OPT_MIN),
            ],
        )
        self._add_tr_row(l_crit, "label.criterion", self.opt_criterion)

        g_ctrl = self._tr_groupbox("optim.group.controls")
        layout.addWidget(g_ctrl)
        l_ctrl = QVBoxLayout(g_ctrl)
        self.btn_opt_run = self._tr_button("button.optimize")
        self.btn_opt_run.setFixedHeight(40)
        self.btn_opt_run.clicked.connect(self.start_optim)
        self.btn_opt_stop = self._tr_button("button.stop")
        self.btn_opt_stop.setEnabled(False)
        self.btn_opt_stop.clicked.connect(self.stop_optim)
        self.btn_opt_save = self._tr_button("button.export_result")
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
        is_optim = (self.tabs.widget(index) is self.tab_optim)
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

    def on_optim_finish(self, success, key, params):
        message = self._message_from_key(key, params)
        self.btn_opt_run.setEnabled(True)
        self.btn_opt_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if key == "err.critical":
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
            self.conc_value.setSuffix(t("unit.per_m3"))
            if self.conc_value.value() <= 0:
                self.conc_value.setValue(1e9)
        else:
            self.conc_value.setDecimals(2)
            self.conc_value.setSuffix(t("unit.g_per_m3"))
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
            self.inv_target_value.setSuffix(t("unit.per_m"))
        else:
            self.inv_target_value.setDecimals(6)
            self.inv_target_value.setSuffix(t("unit.m2_per_g"))

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

    def on_finish(self, success, key, params):
        message = self._message_from_key(key, params)
        self.btn_run.setEnabled(True)
        self.btn_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if key == "err.critical":
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

    def on_inverse_finish(self, success, key, params):
        message = self._message_from_key(key, params)
        self.btn_inv_run.setEnabled(True)
        self.btn_inv_stop.setEnabled(False)
        self.lbl_st.setText(message)

        if success:
            self.pbar.setValue(100)
            QMessageBox.information(self, "Готово", message)
        else:
            self.pbar.setValue(0)
            self.log_msg(f"\n!!! {message} !!!")
            if key == "err.critical":
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
    _set_windows_app_user_model_id()
    QApplication.setOrganizationName(APP_NAME)
    QApplication.setApplicationName(APP_NAME)
    QApplication.setDesktopFileName(APP_NAME)
    app = QApplication(sys.argv)
    i18n.load_from_settings()
    app.setWindowIcon(_application_icon())
    app.setStyle("Fusion")
    w = MainWindow()
    w.show()
    sys.exit(app.exec())
