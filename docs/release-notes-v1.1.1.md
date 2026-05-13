A material-database update for Mie-theory optical attenuation calculations, expanding the built-in refractive-index coverage and documenting model quality more explicitly.

## Downloads

- **macOS (Apple Silicon):** `MieShield-1.1.1-macos-arm64.dmg`
- **Windows (x64):** `MieShield.exe`
- **Linux (x64):** `MieShield-1.1.1-linux-x64.tar.gz`

The Linux archive contains a standalone `mie_shield.dist` folder with the `MieShield` executable and bundled dependencies.

## Highlights

- Expanded the built-in material database from 10 to 26 materials.
- Added `K2CO3`, `KCN`, `AlN`, `MgAl2O4`, `KAlO2`, `Na2CO3`, `ZrO2` (monoclinic baddeleyite), `YSZ`, `ZrC`, `TiC`, `BaO`, `Ba3N2`, `Mg3N2`, `NaAlO2`, `SrO`, and `AlF3`.
- Added bilingual English/Russian material labels for the new entries.
- Documented model quality in the README and expanded the refractive-index handbook with source notes, density policy, wavelength ranges, and reconstruction caveats.
- Added literature-backed primary models where available, including AlN, MgAl2O4, YSZ, and TiC.
- Added DFPT-, handbook-, and Lorentz-reconstruction models for materials without complete open primary `n(lambda), k(lambda)` datasets.
- Kept pure `ZrO2` as monoclinic baddeleyite and exposed cubic yttria-stabilized zirconia separately as `YSZ`.
- Updated the release workflow so GitHub Releases automatically use `docs/release-notes-<tag>.md` for the selected tag.
- Fixed scrollbars in all three tabs: content is now always reachable even with 26 materials listed; the vertical scrollbar no longer overlaps input fields.
- Widened the controls panel (stretch ratio 2:3) and set a minimum tab width of 300 px so the UI remains comfortable at smaller window sizes.

Some new materials are engineering-grade reconstructions rather than measured spectra. See `complex-refractive-indices.md` before using Quality-B or Quality-C entries for publication-grade calculations.
