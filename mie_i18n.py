from __future__ import annotations

from PySide6.QtCore import QObject, QSettings, Signal

from mie_i18n_strings import DEFAULT_LANGUAGE, STRINGS, SUPPORTED_LANGUAGES


class I18n(QObject):
    language_changed = Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self._language = DEFAULT_LANGUAGE

    @property
    def language(self) -> str:
        return self._language

    def load_from_settings(self) -> None:
        lang = str(QSettings().value("ui/language", DEFAULT_LANGUAGE))
        self._language = lang if lang in SUPPORTED_LANGUAGES else DEFAULT_LANGUAGE

    def set_language(self, language: str) -> None:
        lang = language if language in SUPPORTED_LANGUAGES else DEFAULT_LANGUAGE
        if lang == self._language:
            return
        self._language = lang
        QSettings().setValue("ui/language", lang)
        self.language_changed.emit(lang)

    def t(self, key: str, **kwargs: object) -> str:
        text = STRINGS.get(self._language, {}).get(key)
        if text is None:
            text = STRINGS[DEFAULT_LANGUAGE].get(key, key)
        if kwargs:
            return text.format(**kwargs)
        return text


i18n = I18n()


def t(key: str, **kwargs: object) -> str:
    return i18n.t(key, **kwargs)
