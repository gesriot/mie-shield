import string

from PySide6.QtWidgets import QApplication

from mie_i18n import i18n
from mie_i18n_strings import STRINGS
from mie_shield import MainWindow


def test_i18n_catalog_keys_match():
    assert set(STRINGS["en"]) == set(STRINGS["ru"])


def test_i18n_catalog_placeholders_match():
    formatter = string.Formatter()
    for key in STRINGS["en"]:
        en_fields = {field for _, field, _, _ in formatter.parse(STRINGS["en"][key]) if field}
        ru_fields = {field for _, field, _, _ in formatter.parse(STRINGS["ru"][key]) if field}
        assert en_fields == ru_fields, key


def test_main_window_live_switches_tab_labels():
    app = QApplication.instance() or QApplication([])
    original_language = i18n.language
    i18n.set_language("en")
    window = MainWindow()
    try:
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_forward)) == "Forward problem"
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_inverse)) == "Inverse problem"
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_optim)) == "Optimization"
        assert window.dist_type_combo.itemText(0) == "Monodisperse"
        assert window.btn_run.text() == "RUN"
        assert window.mats["C"][0].text() == "C - Carbon (Soot)"
        assert window.lbl_st.text() == "System ready."

        i18n.set_language("ru")
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_forward)) == "Прямая задача"
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_inverse)) == "Обратная задача"
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_optim)) == "Оптимизация"
        assert window.dist_type_combo.itemText(0) == "Монодисперс"
        assert window.btn_run.text() == "ЗАПУСК"
        assert window.mats["C"][0].text() == "C - Углерод (Soot)"
        assert window.lbl_st.text() == "Система готова."

        i18n.set_language("en")
        assert window.tabs.tabText(window.tabs.indexOf(window.tab_forward)) == "Forward problem"
        assert window.lbl_st.text() == "System ready."
    finally:
        i18n.set_language(original_language)
        window.close()
        app.processEvents()
