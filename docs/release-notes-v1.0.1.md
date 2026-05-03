A desktop GUI update for Mie-theory optical attenuation calculations, focused on bilingual UI support and release reliability.

## Downloads

- **macOS (Apple Silicon):** `MieShield-1.0.1-macos-arm64.dmg`
- **Windows (x64):** `MieShield.exe`
- **Linux (x64):** `MieShield-1.0.1-linux-x64.tar.gz`

The Linux archive contains a standalone `mie_shield.dist` folder with the `MieShield` executable and bundled dependencies.

## Highlights

- Added a bilingual desktop UI: English is the default, and Russian is available from the `Language` menu.
- Language switching is live and persists across launches via `QSettings`.
- Translation strings are stored in Python modules, so Nuitka onefile/standalone builds do not need external localization files.
- Localized tabs, menus, controls, material display names, validation dialogs, runtime statuses, and standard Yes/No confirmation buttons.
- Kept calculation logs and exported reports as stable English technical output.
- Made core validation errors language-neutral with stable error codes and UI-side rendering.
- Added i18n regression coverage for catalog parity, placeholder parity, live tab/control switching, and localized material names.
- Hardened release artifact uploads so missing macOS, Windows, or Linux assets fail the workflow instead of warning.
