# mie-shield

Mie scattering calculator for extinction, MEC, and transmission in airborne particle mixtures.

`mie-shield` is a desktop GUI for estimating optical attenuation by dispersed particles using Mie theory. It calculates extinction cross sections, attenuation coefficients, optical depth, transmittance, and mass extinction coefficient (MEC) for monodisperse and polydisperse mixtures.

The application is intended for engineering and experimental analysis of aerodisperse particle media: soot, metals, salts, oxides, and mixed particles with user-defined number fractions.

## Features

- Forward Mie calculation for monodisperse particles.
- Forward Mie calculation for log-normal and custom log-normal-like particle-size distributions.
- Mixtures of several particle materials by number fraction.
- Mass or number concentration input.
- Configurable measurement path length `L`.
- Spectral calculations over a wavelength range.
- Inverse diameter search from a target `MEC`, `alpha`, `tau`, `T_eff`, or `AVG T`.
- Optimization of a custom particle-size distribution window for maximum `MEC * L`.
- Export of forward, inverse, and optimization results to text files.

## Materials

The built-in material database currently includes:

| Code | Material |
| --- | --- |
| `C` | carbon soot |
| `Mg` | magnesium |
| `MgCl2` | magnesium chloride |
| `ZnCl2` | zinc chloride |
| `MgF2` | magnesium fluoride |
| `Al4C3` | aluminum carbide |
| `Al` | aluminum |
| `MgO` | magnesium oxide |
| `Al2O3` | alumina / corundum |
| `CuZn` | Cu70/Zn30 brass |

Optical constants are documented in [complex-refractive-indices.md](complex-refractive-indices.md). Some materials use literature-based models, while others use approximate reconstructed models where complete measured `n(lambda), k(lambda)` tables are unavailable. Treat approximate materials accordingly.

## Installation

The project uses [`uv`](https://docs.astral.sh/uv/) and includes a lock file.

```bash
uv sync --locked --group dev
```

## Running

```bash
uv run --locked python mie_shield.py
```

The main window has three calculation tabs:

- **Forward problem** (`Прямая задача`)
- **Inverse problem** (`Обратная задача`)
- **Optimization** (`Оптимизация`)

## Installation

The project uses [`uv`](https://docs.astral.sh/uv/) and includes a lock file.

```bash
uv sync --locked --group dev
```

## Running

```bash
uv run --locked python mie_shield.py
```

The main window has three calculation tabs:

- **Forward problem** (`Прямая задача`)
- **Inverse problem** (`Обратная задача`)
- **Optimization** (`Оптимизация`)

## Packaging

Build tools are intentionally kept out of the project runtime dependencies.
The application depends on PySide6, NumPy, SciPy, Matplotlib, and PyMieScatt;
packaging additionally uses Nuitka, but only inside a local build environment
created by the scripts under `scripts/`.

Generated artifacts are written to `dist/` and should not be committed. Attach
DMG/ZIP files to GitHub Releases instead.

### macOS DMG

Requirements:

- macOS on the target architecture.
- `uv`.
- Python 3.14 available to `uv`.
- Xcode Command Line Tools (`clang`, `codesign`, `hdiutil`).

Build:

```bash
./scripts/package-macos.sh
```

The script creates `.build-venv/macos`, installs runtime dependencies plus
build-only Nuitka there, builds a signed ad-hoc `.app`, and packages it into:

```text
dist/MieShield-<version>-macos-arm64.dmg
```

If `icon.png` is present in the repository root, the script crops it to the
visible rounded-square frame, creates transparent macOS icon corners, and
embeds the generated `.icns` in the app bundle.

Useful overrides:

```bash
PYTHON_VERSION=3.14 ARCH=arm64 PRODUCT_NAME=MieShield ./scripts/package-macos.sh
```

The release build uses Nuitka with link-time optimization and optimized Python
runtime flags:

```text
--lto=yes --python-flag=-O --python-flag=no_docstrings
```

### Windows ZIP

Windows builds should be produced on Windows, not cross-compiled from macOS.
The first Windows packaging target is a standalone Nuitka folder zipped for
release distribution.

Requirements:

- Windows.
- `uv`.
- Python 3.14 available to `uv`.
- Microsoft C++ Build Tools.

Build from PowerShell:

```powershell
.\scripts\package-windows.ps1
```

The script creates `.build-venv/windows`, installs build-only Nuitka there,
builds a GUI executable without a console window, and writes:

```text
dist/MieShield-<version>-windows-x64.zip
```

An installer can be added later with Inno Setup, WiX, or NSIS while keeping the
same Nuitka build output as input.

## Core Quantities

For monodisperse particles, the main quantities are:

```text
Cext = Qext * pi * D^2 / 4
alpha = N * Cext
tau = alpha * L
T = exp(-tau)
MEC = alpha / rho_mass
```

where:

- `D` is particle diameter.
- `Qext` is the Mie extinction efficiency.
- `Cext` is extinction cross section.
- `N` is number concentration.
- `rho_mass` is mass concentration in `g/m^3`.
- `L` is the measurement path length in meters.
- `tau` is optical depth.
- `T` is transmittance.
- `MEC` is mass extinction coefficient in `m^2/g`.

For wavelength ranges, the app reports both:

- `AVG T = mean(exp(-alpha_i * L))`
- `T_eff = exp(-mean(alpha_i * L))`

These are not the same in general. Use `AVG T` when the experimental observable is an arithmetic average of spectral transmittance. Use `T_eff` when the observable is represented by an effective optical depth.

## Limitations

- Particles are modeled as spheres using Mie theory.
- Mixture entries are treated as number-fraction weights and normalized internally.
- Some refractive-index models are approximate; see the optical-constant documentation before using results as publication-grade material constants.
- The optimization tab maximizes the selected model metric. It does not infer a unique physical particle-size distribution from experimental data.

## License

MIT – see [LICENSE](LICENSE).
