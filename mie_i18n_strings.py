from __future__ import annotations

DEFAULT_LANGUAGE = "en"
SUPPORTED_LANGUAGES = ("en", "ru")

LANGUAGE_NAMES = {
    "en": "English",
    "ru": "Русский",
}

STRINGS = {
    "en": {
        "app.title": "Mie Extinction Calculator (Alpha + MEC)",
        "menu.language": "Language",
        "status.ready": "System ready.",
        "tab.forward": "Forward problem",
        "tab.inverse": "Inverse problem",
        "tab.optim": "Optimization",
        "material.C": "Carbon (Soot)",
        "material.Mg": "Magnesium (Metal)",
        "material.MgCl2": "Magnesium chloride",
        "material.ZnCl2": "Zinc chloride",
        "material.MgF2": "Magnesium fluoride",
        "material.Al4C3": "Aluminum carbide",
        "material.Al": "Aluminum (Metal)",
        "material.MgO": "Magnesium oxide",
        "material.Al2O3": "Aluminum oxide (corundum)",
        "material.CuZn": "Brass Cu70/Zn30 (Brass)",
        "err.conc_invalid_mode": "Unknown concentration mode: {conc_mode!r}.",
        "err.conc_non_positive": "Concentration must be > 0.",
        "err.avg_mass_invalid": "Average particle mass is invalid (<= 0 or NaN).",
        "err.number_conc_invalid": "Number concentration is invalid (<= 0 or NaN).",
        "err.mass_conc_positive": "Mass concentration must be > 0.",
        "err.transmittance_out_of_range": "Transmittance T must be in the range (0, 1].",
    },
    "ru": {
        "app.title": "Mie Extinction Calculator (Alpha + MEC)",
        "menu.language": "Язык",
        "status.ready": "Система готова.",
        "tab.forward": "Прямая задача",
        "tab.inverse": "Обратная задача",
        "tab.optim": "Оптимизация",
        "material.C": "Углерод (Soot)",
        "material.Mg": "Магний (Metal)",
        "material.MgCl2": "Хлорид магния",
        "material.ZnCl2": "Хлорид цинка",
        "material.MgF2": "Фторид магния",
        "material.Al4C3": "Карбид алюминия",
        "material.Al": "Алюминий (Metal)",
        "material.MgO": "Оксид магния",
        "material.Al2O3": "Оксид алюминия (corundum)",
        "material.CuZn": "Латунь Cu70/Zn30 (Brass)",
        "err.conc_invalid_mode": "Неизвестный режим концентрации: {conc_mode!r}.",
        "err.conc_non_positive": "Концентрация должна быть > 0.",
        "err.avg_mass_invalid": "Средняя масса частицы некорректна (<= 0 или NaN).",
        "err.number_conc_invalid": "Числовая концентрация некорректна (<= 0 или NaN).",
        "err.mass_conc_positive": "Массовая концентрация должна быть > 0.",
        "err.transmittance_out_of_range": "Пропускание T должно быть в диапазоне (0, 1].",
    },
}
