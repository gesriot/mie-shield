import math
from types import SimpleNamespace

import numpy as np
import pytest

import mie_core as mc
import mie_shield as ms

# Patch dependencies in the module where the tested function resolves them:
# core helpers use mc.*, shield workers use ms.* globals.


def run_worker(worker):
    logs = []
    progress = []
    results = []
    finished = []

    worker.log_signal.connect(logs.append)
    worker.progress_signal.connect(progress.append)
    worker.result_signal.connect(results.append)
    worker.finished_signal.connect(lambda success, message: finished.append((success, message)))

    worker.run()

    assert finished, "worker did not emit finished_signal"
    return SimpleNamespace(
        logs=logs,
        progress=progress,
        results=results,
        finished=finished[-1],
    )


def fake_mieq(q_ext):
    def _fake_mieq(_m, _wavelength_nm, _diameter_nm, asDict=False):
        if asDict:
            return {"Qext": q_ext}
        return (q_ext, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    return _fake_mieq


def mono_params(**overrides):
    params = {
        "fractions": {"C": 1.0},
        "wl_range": (0.5, 0.6, 0.1),
        "path_length_m": 2.7,
        "conc_mode": ms.CONC_MASS,
        "conc_value": 1.2,
        "monodisperse": True,
        "D_um": 1.0,
    }
    params.update(overrides)
    return params


def distribution_params(dist_type, **overrides):
    params = {
        "fractions": {"C": 0.7, "ZnCl2": 0.3},
        "wl_range": (0.5, 0.6, 0.1),
        "path_length_m": 2.7,
        "conc_mode": ms.CONC_MASS,
        "conc_value": 1.2,
        "monodisperse": False,
        "d_range": (0.1, 1.0),
        "points_d": 9,
        "dist_type": dist_type,
        "d_dist": (0.4, 1.5),
        "custom_A": 1.0,
        "custom_mu": -0.5,
        "custom_sigma": 0.4,
    }
    params.update(overrides)
    return params


def expected_mono_values(q_ext=2.0, rho=1800.0, d_um=1.0, mass_conc_g=1.2, path_length_m=2.7):
    d_m = d_um * 1e-6
    d_nm = d_um * 1000.0
    avg_mass_kg = rho * (math.pi / 6.0) * d_m**3
    num_conc = mass_conc_g * 1e-3 / avg_mass_kg
    cext_um2 = q_ext * (math.pi * d_nm**2 / 4.0) / 1e6
    alpha = num_conc * cext_um2 * 1e-12
    tau = alpha * path_length_m
    transmittance = math.exp(-tau)
    mec = alpha / mass_conc_g
    return SimpleNamespace(
        avg_mass_kg=avg_mass_kg,
        num_conc=num_conc,
        mass_conc_g=mass_conc_g,
        cext_um2=cext_um2,
        alpha=alpha,
        tau=tau,
        transmittance=transmittance,
        mec=mec,
    )


def test_make_wavelengths_uses_endpoints_and_actual_step():
    wavelengths = ms.make_wavelengths(0.4, 0.78, 0.1)

    assert wavelengths.tolist() == pytest.approx([0.4, 0.495, 0.59, 0.685, 0.78])
    assert ms.make_wavelengths(1.0, 0.5, 0.1).size == 0
    assert ms.make_wavelengths(0.5, 1.0, 0.0).size == 0


def test_refractive_indices_are_finite_for_all_materials():
    material_codes = [code for code, _name, _rho in ms.MATERIALS_DB]

    for material in material_codes:
        for lam_um in (0.4, 0.59, 5.0, 10.0, 20.0, 50.0):
            m = ms.get_ri(material, lam_um)
            assert np.isfinite(m.real), (material, lam_um, m)
            assert np.isfinite(m.imag), (material, lam_um, m)
            assert m.real > 0.0, (material, lam_um, m)
            assert m.imag >= 0.0, (material, lam_um, m)


def test_zincl2_visible_index_and_density_are_in_database():
    density_by_code = {code: rho for code, _name, rho in ms.MATERIALS_DB}
    m = ms.get_ri("ZnCl2", 0.589)

    assert density_by_code["ZnCl2"] == pytest.approx(2907.0)
    assert m.real == pytest.approx(math.sqrt(2.87), rel=1e-2)
    assert m.imag < 1e-3


def test_qext_to_cext_um2_uses_geometric_cross_section():
    assert ms.qext_to_cext_um2(2.0, 1000.0) == pytest.approx(math.pi / 2.0)


def test_safe_mie_qext_matches_direct_mieq():
    lam_nm = 500.0
    D_nm = 1000.0
    m = ms.get_ri("C", lam_nm / 1000.0)
    expected = mc.MieQ(m, lam_nm, D_nm, asDict=False)[0]

    q_ext, ok = ms.safe_mie_qext(m, lam_nm, D_nm)

    assert ok is True
    assert q_ext == pytest.approx(expected)


def test_safe_mie_qext_returns_failure_on_exception(monkeypatch):
    def failing_mieq(_m, _wavelength_nm, _diameter_nm, asDict=False):
        raise RuntimeError("forced Mie failure")

    monkeypatch.setattr(mc, "MieQ", failing_mieq)

    assert ms.safe_mie_qext(1.5 + 0.0j, 500.0, 1000.0) == (0.0, False)


@pytest.mark.parametrize(
    ("q_ext", "expected"),
    [
        (np.nan, (0.0, False)),
        (-1e-8, (0.0, False)),
        (-5e-10, (0.0, True)),
    ],
)
def test_safe_mie_qext_handles_invalid_and_small_negative_values(monkeypatch, q_ext, expected):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(q_ext))

    assert ms.safe_mie_qext(1.5 + 0.0j, 500.0, 1000.0) == expected


def test_conc_mode_constants_are_machine_strings():
    # ConcMode values are decoupled from UI labels (Russian "Массовая"/"Числовая"
    # live only in the QComboBox addItem() text; userData holds these constants).
    assert ms.CONC_MASS == "mass"
    assert ms.CONC_NUMBER == "number"


def test_resolve_concentration_rejects_unknown_mode():
    # Guards against silent fall-through when a non-canonical string slips in
    # (e.g. an old UI label leaking into a worker dict). Without this check
    # an unknown mode was silently treated as mass concentration.
    with pytest.raises(ValueError, match="Неизвестный режим концентрации"):
        ms.resolve_concentration("Числовая", 5e8, 1e-15)


def test_resolve_concentration_round_trip_mass_and_number():
    avg_mass_kg = 1.0e-15  # arbitrary positive
    mass_conc_g = 1.5

    num_from_mass, mass_back = ms.resolve_concentration(ms.CONC_MASS, mass_conc_g, avg_mass_kg)
    expected_num = (mass_conc_g * 1e-3) / avg_mass_kg
    assert num_from_mass == pytest.approx(expected_num)
    assert mass_back == pytest.approx(mass_conc_g)

    num_back, mass_from_num = ms.resolve_concentration(ms.CONC_NUMBER, num_from_mass, avg_mass_kg)
    assert num_back == pytest.approx(expected_num)
    assert mass_from_num == pytest.approx(mass_conc_g)


def test_compute_qext_avg_single_material(monkeypatch):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(1.25))

    qext_avg, lost_weight = ms.compute_qext_avg({"C": 1.0}, 0.5, 1.0)

    assert qext_avg == pytest.approx(1.25)
    assert lost_weight == pytest.approx(0.0)


def test_compute_qext_avg_weights_material_qext(monkeypatch):
    def fake_get_ri(material, _lam_um):
        return material

    def material_mieq(material, _wavelength_nm, _diameter_nm, asDict=False):
        return ({"C": 1.0, "ZnCl2": 3.0}[material], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    monkeypatch.setattr(mc, "get_ri", fake_get_ri)
    monkeypatch.setattr(mc, "MieQ", material_mieq)

    qext_avg, lost_weight = ms.compute_qext_avg({"C": 0.25, "ZnCl2": 0.75}, 0.5, 1.0)

    assert qext_avg == pytest.approx(2.5)
    assert lost_weight == pytest.approx(0.0)


def test_compute_qext_avg_accumulates_lost_weight(monkeypatch):
    def failing_mieq(_m, _wavelength_nm, _diameter_nm, asDict=False):
        raise RuntimeError("forced Mie failure")

    monkeypatch.setattr(mc, "MieQ", failing_mieq)

    qext_avg, lost_weight = ms.compute_qext_avg({"C": 0.4, "ZnCl2": 0.6}, 0.5, 1.0)

    assert qext_avg == pytest.approx(0.0)
    assert lost_weight == pytest.approx(1.0)


def test_compute_mec_for_d_uses_qext_average(monkeypatch):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(2.0))
    fractions = {"C": 0.25, "ZnCl2": 0.75}
    rho_avg = 0.25 * 1800.0 + 0.75 * 2907.0

    mec = ms.compute_mec_for_d(1.0, fractions, rho_avg, 0.5)

    expected = (3.0 / 2.0) * (1e-3 / rho_avg) * 2.0 / 1e-6
    assert mec == pytest.approx(expected)


def test_compute_mec_for_d_returns_nan_when_too_much_mie_weight_is_lost(monkeypatch):
    def failing_mieq(_m, _wavelength_nm, _diameter_nm, asDict=False):
        raise RuntimeError("forced Mie failure")

    monkeypatch.setattr(mc, "MieQ", failing_mieq)

    mec = ms.compute_mec_for_d(1.0, {"C": 1.0}, 1800.0, 0.5)

    assert np.isnan(mec)


def test_nan_safe_bisect_handles_regular_and_nan_midpoints():
    assert ms.nan_safe_bisect(lambda x: x - 2.0, 1.0, 3.0) == pytest.approx(2.0, abs=1e-6)

    def f_with_nan_gap(x):
        if 1.45 < x < 1.55:
            return np.nan
        return x - 1.25

    assert ms.nan_safe_bisect(f_with_nan_gap, 1.0, 2.0) == pytest.approx(1.25, abs=1e-5)


def test_calculation_worker_monodisperse_mass_concentration_uses_path_length(monkeypatch):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(2.0))

    out = run_worker(ms.CalculationWorker(mono_params()))

    assert out.finished == (True, "Расчет успешно завершен.")
    result = out.results[-1]
    expected = expected_mono_values()
    assert result["avg_mass_kg"] == pytest.approx(expected.avg_mass_kg)
    assert result["num_conc"] == pytest.approx(expected.num_conc)
    assert result["mass_conc_g"] == pytest.approx(expected.mass_conc_g)

    for row in result["data"]:
        assert row["cext_um2"] == pytest.approx(expected.cext_um2)
        assert row["alpha_1m"] == pytest.approx(expected.alpha)
        assert row["tau"] == pytest.approx(expected.tau)
        assert row["transmittance"] == pytest.approx(expected.transmittance)
        assert row["mec_m2g"] == pytest.approx(expected.mec)
        assert row["parts"]["C"] == pytest.approx(expected.cext_um2)

    assert any("T_eff=exp(-AVG tau)" in line for line in out.logs)


def test_calculation_worker_monodisperse_number_concentration(monkeypatch):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(2.0))
    number_conc = 5.0e8
    expected = expected_mono_values(mass_conc_g=1.0)
    expected_mass_conc_g = number_conc * expected.avg_mass_kg * 1000.0
    expected_alpha = number_conc * expected.cext_um2 * 1e-12

    out = run_worker(
        ms.CalculationWorker(
            mono_params(conc_mode=ms.CONC_NUMBER, conc_value=number_conc)
        )
    )

    assert out.finished[0] is True
    result = out.results[-1]
    assert result["num_conc"] == pytest.approx(number_conc)
    assert result["mass_conc_g"] == pytest.approx(expected_mass_conc_g)
    assert result["data"][0]["alpha_1m"] == pytest.approx(expected_alpha)
    assert result["data"][0]["mec_m2g"] == pytest.approx(expected_alpha / expected_mass_conc_g)


@pytest.mark.parametrize(
    ("dist_type", "expected_log"),
    [
        ("lognormal", "Truncated Log-Normal"),
        ("custom", "Кастомное распределение"),
    ],
)
def test_calculation_worker_polydisperse_distribution_branches(monkeypatch, dist_type, expected_log):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(2.0))

    out = run_worker(ms.CalculationWorker(distribution_params(dist_type)))

    assert out.finished[0] is True
    assert any(expected_log in line for line in out.logs)
    result = out.results[-1]
    assert result["mass_conc_g"] == pytest.approx(1.2)
    assert result["avg_mass_kg"] > 0.0
    assert len(result["data"]) == 2
    for row in result["data"]:
        assert row["alpha_1m"] > 0.0
        assert row["tau"] == pytest.approx(row["alpha_1m"] * 2.7)
        assert row["transmittance"] == pytest.approx(math.exp(-row["tau"]))
        assert row["mec_m2g"] == pytest.approx(row["alpha_1m"] / 1.2)


def test_calculation_worker_rejects_invalid_path_length(monkeypatch):
    monkeypatch.setattr(mc, "MieQ", fake_mieq(2.0))

    out = run_worker(ms.CalculationWorker(mono_params(path_length_m=0.0)))

    assert out.finished[0] is False
    assert "Длина трассы L должна быть > 0" in out.finished[1]


def inverse_params(input_mode, target_value, wl_mode="single", **overrides):
    params = {
        "fractions": {"C": 1.0},
        "wl_mode": wl_mode,
        "lambda_um": 0.5,
        "wl_range": (0.4, 0.6, 0.1),
        "input_mode": input_mode,
        "target_value": target_value,
        "mass_conc_g": 1.2,
        "path_length_m": 2.7,
        "D_min_um": 0.2,
        "D_max_um": 5.0,
        "N_scan": 80,
    }
    params.update(overrides)
    return params


@pytest.mark.parametrize(
    ("input_mode", "target_value", "wl_mode", "expected_label", "expected_metric"),
    [
        ("mec", 2.0, "single", "MEC", 2.0),
        ("mec", 2.0, "range", "AVG MEC", 2.0),
        ("alpha", 2.4, "single", "MEC", 2.0),
        ("tau", 6.48, "single", "MEC", 2.0),
        ("effective_transmittance", math.exp(-6.48), "range", "T_eff", math.exp(-6.48)),
        ("transmittance", math.exp(-6.48), "range", "AVG T", math.exp(-6.48)),
    ],
)
def test_inverse_worker_solves_input_modes(monkeypatch, input_mode, target_value, wl_mode, expected_label, expected_metric):
    def fake_compute_mec_for_d(d_um, _fractions, _rho_avg, _lam_um):
        return 2.0 / d_um

    monkeypatch.setattr(ms, "compute_mec_for_d", fake_compute_mec_for_d)

    out = run_worker(ms.InverseWorker(inverse_params(input_mode, target_value, wl_mode)))

    assert out.finished[0] is True
    result = out.results[-1]
    assert result["metric_label"] == expected_label
    assert result["solutions"], out.logs
    solution = result["solutions"][0]
    assert solution["D_um"] == pytest.approx(1.0, abs=2e-4)
    assert solution["mec_check"] == pytest.approx(expected_metric, rel=2e-4)

    if input_mode in ("alpha", "tau", "transmittance", "effective_transmittance"):
        assert solution["N"] is not None
        assert solution["avg_mec_check"] == pytest.approx(2.0, rel=2e-4)
        assert solution["alpha_check"] == pytest.approx(2.4, rel=2e-4)
        assert solution["tau_check"] == pytest.approx(6.48, rel=2e-4)
        assert solution["transmittance_check"] == pytest.approx(math.exp(-6.48), rel=2e-4)

    if input_mode in ("transmittance", "effective_transmittance"):
        assert any("AVG MEC (справка)" in line for line in out.logs)

    if input_mode == "transmittance" and wl_mode == "range":
        assert any("exp(-AVG τ) (справка)" in line for line in out.logs)
    if input_mode == "effective_transmittance" and wl_mode == "range":
        assert any("exp(-AVG τ) (проверка)" in line for line in out.logs)


def test_inverse_worker_reports_no_solution(monkeypatch):
    monkeypatch.setattr(ms, "compute_mec_for_d", lambda d_um, *_args: 2.0 / d_um)

    out = run_worker(ms.InverseWorker(inverse_params("mec", 100.0)))

    assert out.finished[0] is True
    assert out.results[-1]["solutions"] == []
    assert any("РЕШЕНИЙ НЕ НАЙДЕНО" in line for line in out.logs)


def test_inverse_worker_rejects_invalid_transmittance(monkeypatch):
    monkeypatch.setattr(ms, "compute_mec_for_d", lambda d_um, *_args: 2.0 / d_um)

    out = run_worker(ms.InverseWorker(inverse_params("effective_transmittance", 1.1)))

    assert out.finished[0] is False
    assert "Пропускание T должно быть" in out.finished[1]


def optimization_params(mode, criterion):
    return {
        "fractions": {"C": 1.0},
        "wl_range": (0.4, 0.5, 0.1),
        "path_length_m": 2.7,
        "D_scan_min": 0.1,
        "D_scan_max": 1.0,
        "N_D_scan": 5,
        "N_D_points": 7,
        "N_window_grid": 4,
        "criterion": criterion,
        "mode": mode,
        "mu_fixed": 0.0,
        "sigma_fixed": 0.5,
        "mu_range": (-0.2, 0.2),
        "sigma_range": (0.2, 0.8),
    }


@pytest.mark.parametrize(
    ("mode", "criterion"),
    [
        ("window_only", "mean"),
        ("window_only", "min"),
        ("full", "mean"),
        ("full", "min"),
    ],
)
def test_optimization_worker_modes_and_criteria(monkeypatch, mode, criterion):
    def fake_compute_mec_for_d(d_um, _fractions, _rho_avg, lam_um):
        return (1.0 + 0.1 * lam_um) / (d_um + 0.1)

    def fake_differential_evolution(objective, bounds, **_kwargs):
        x = np.array([(lo + hi) / 2.0 for lo, hi in bounds], dtype=float)
        value = objective(x)
        assert np.isfinite(value)
        return SimpleNamespace(x=x)

    monkeypatch.setattr(ms, "compute_mec_for_d", fake_compute_mec_for_d)
    monkeypatch.setattr(ms.scipy.optimize, "differential_evolution", fake_differential_evolution)

    out = run_worker(ms.OptimizationWorker(optimization_params(mode, criterion)))

    assert out.finished[0] is True
    heatmap = next(result for result in out.results if result["type"] == "heatmap")
    window_heatmap = next(result for result in out.results if result["type"] == "window_heatmap")
    result = next(result for result in out.results if result["type"] == "result")

    assert heatmap["metric_label"] == "MEC·L"
    assert np.array(heatmap["path_score_map"]) == pytest.approx(np.array(heatmap["mec_map"]) * 2.7)
    assert window_heatmap["metric_label"] == "MEC·L"
    assert window_heatmap["criterion"] == criterion
    assert result["criterion"] == criterion
    assert np.array(result["mec_l_spectrum"]) == pytest.approx(np.array(result["mec_spectrum"]) * 2.7)
    assert result["mec_l_mean"] == pytest.approx(result["mec_mean"] * 2.7)
    assert result["mec_l_min"] == pytest.approx(result["mec_min"] * 2.7)


def test_optimization_worker_rejects_empty_spectrum(monkeypatch):
    monkeypatch.setattr(ms, "compute_mec_for_d", lambda d_um, *_args: 1.0 / d_um)

    out = run_worker(
        ms.OptimizationWorker(
            optimization_params("window_only", "mean") | {"wl_range": (0.5, 0.4, 0.1)}
        )
    )

    assert out.finished == (False, "Пустой спектр длин волн.")
