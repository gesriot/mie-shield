A desktop GUI for Mie-theory calculations of optical attenuation by aerodisperse particle mixtures (extinction, MEC, optical depth, transmittance).

## Downloads

- **macOS (Apple Silicon):** `MieShield-1.0.0-macos-arm64.dmg`
- **Windows (x64):** `MieShield.exe`
- **Linux (x64):** `MieShield-1.0.0-linux-x64.tar.gz`

The Linux archive contains a standalone `mie_shield.dist` folder with the `MieShield` executable and bundled dependencies.

## Highlights

- Forward Mie calculation for monodisperse particles and log-normal / custom distributions.
- Mixtures of multiple materials by number fraction with mass or number concentration input.
- Spectral calculations over a wavelength range, with both `AVG T` and `T_eff` reported.
- Inverse diameter search from a target `MEC`, `alpha`, `tau`, `T_eff`, or `AVG T`.
- Optimization of a custom particle-size distribution window for maximum `MEC * L`.
- Export of forward, inverse, and optimization results to text files.
- Built-in material database: `C`, `Mg`, `MgCl2`, `ZnCl2`, `MgF2`, `Al4C3`, `Al`, `MgO`, `Al2O3`, `CuZn` - see [complex-refractive-indices.md](https://github.com/gesriot/mie-shield/blob/main/complex-refractive-indices.md) for sources and which models are approximate.
