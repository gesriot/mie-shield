# Complex refractive indices for Mie scattering: a 10-material reference

This code combines literature-based optical-constant models with approximate or reconstructed models for materials where complete measured n(λ),k(λ) datasets are unavailable. Several entries reproduce published coefficients or tabulated values directly, while MgCl₂, ZnCl₂, Mg, and especially Al₄C₃ should be treated as engineering approximations with explicit caveats. This guide documents each implemented model, flags discrepancies, and records wavelength ranges and limitations relevant for Mie scattering calculations. The Sellmeier C-values throughout require careful attention: for MgF₂, Al₂O₃, and Al₄C₃ they represent resonance *wavelengths* in µm (not λ²), a distinction that produces ~4% refractive index errors if mishandled. For ZnCl₂, where no measured n(λ),k(λ) tables exist, the IR dispersion is represented by a four-oscillator Lorentz reconstruction constrained by first-principles DFPT dielectric tensors (Materials Project mp-22909) and by FTIR/Raman phonon-mode positions from Angell & Wong (1970) and Janz & James (1974).

---

## 1. Carbon soot tracks Chang & Charalampopoulos (1990) with minor visible-band deviations

**Source confirmed:** H. Chang and T. T. Charalampopoulos, "Determination of the wavelength dependence of refractive indices of flame soot," *Proc. R. Soc. London A* **430**, 577–591 (1990).

**Mathematical model – cubic polynomial in ln(λ):**

| Coefficient | n polynomial | k polynomial |
|:-----------:|:------------:|:------------:|
| a₀ | 1.811 | 0.5821 |
| a₁ | 0.1263 | 0.1213 |
| a₂ | 0.0270 | 0.2309 |
| a₃ | 0.0417 | −0.01 |

The formula is n(λ) = a₀ + a₁·ln λ + a₂·(ln λ)² + a₃·(ln λ)³ (and similarly for k), where **λ is in micrometers**. These coefficients are widely cited in combustion modeling literature referencing Chang & Charalampopoulos (1990). The Thermopedia entry by L. A. Dombrovsky confirms the polynomial-in-ln(λ) form with a valid range of **0.4–30 µm**.

A spot-check at λ = 0.55 µm yields n ≈ 1.74, k ≈ 0.59 from the polynomial, while several secondary sources cite the paper's tabulated measurement as **m ≈ 1.77 + 0.63i** at 550 nm. This ~2% discrepancy in n and ~6% in k reflects the polynomial being a least-squares fit rather than exact interpolation of the measured data points. The polynomial performs better in the near-IR through mid-IR where soot radiative properties matter most for heat transfer calculations.

**This is not the Dalzell & Sarofim (1969) formula.** Dalzell and Sarofim used a Drude-Lorentz dispersion model fitted to reflectivity measurements from compressed soot pellets (which may have contained ~⅓ air voids, potentially underestimating n and k by ~20%). Chang & Charalampopoulos used in-situ extinction and scattering measurements on propane flame soot (φ = 1.8, fuel-rich premixed propane-oxygen flame), combined with Kramers-Kronig analysis.

**Key caveats:** These constants are specific to propane soot. Bond & Bergstrom (2006) recommended the higher value **m = 1.95 + 0.79i** for "void-free" black carbon, now preferred in many climate models. Temperature dependence is minimal (~1% per 100 K per Stagg & Charalampopoulos, 1993). Soot maturity, fuel type, and aggregate morphology all affect real optical properties significantly.

---

## 2. Aluminum parameters from Rakić (1998) match to every decimal place

**Source confirmed:** A. D. Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, "Optical properties of metallic films for vertical-cavity optoelectronic devices," *Appl. Opt.* **37**(22), 5271–5283 (1998).

The Lorentz-Drude dielectric function is:

$$\varepsilon(\omega) = 1 - \frac{f_0\,\omega_p^2}{\omega(\omega + i\Gamma_0)} + \sum_{j=1}^{4} \frac{f_j\,\omega_p^2}{\omega_j^2 - \omega^2 - i\omega\Gamma_j}$$

All **15 parameters match the published values exactly**:

| Parameter | Code | Rakić 1998 | Unit |
|:---------:|:----:|:----------:|:----:|
| ωₚ | 14.98 | 14.98 | eV |
| f₀ | 0.523 | 0.523 | – |
| Γ₀ | 0.047 | 0.047 | eV |
| f₁, Γ₁, ω₁ | 0.227, 0.333, 0.162 | 0.227, 0.333, 0.162 | –, eV, eV |
| f₂, Γ₂, ω₂ | 0.050, 0.312, 1.544 | 0.050, 0.312, 1.544 | –, eV, eV |
| f₃, Γ₃, ω₃ | 0.166, 1.351, 1.808 | 0.166, 1.351, 1.808 | –, eV, eV |
| f₄, Γ₄, ω₄ | 0.030, 3.382, 3.473 | 0.030, 3.382, 3.473 | –, eV, eV |

**Valid fitting range: 0.1–6 eV (0.207–12.4 µm).** The model can be evaluated outside this range but reliability degrades. The refractiveindex.info database lists an extended calculation range of 0.0620–248 µm. The Rakić model describes bulk, clean metallic Al at room temperature. Real aluminum surfaces carry a **2–4 nm native Al₂O₃ oxide layer** that affects measured reflectance. Caution is warranted near Re(ε) ≈ 0, where the Drude-type fit can be qualitatively incorrect.

---

## 3. Magnesium uses a coarse approximation to broadband optical data

**Source confirmed:** H.-J. Hagemann, W. Gudat, and C. Kunz, "Optical constants from the far infrared to the x-ray region: Mg, Al, Cu, Ag, Au, Bi, C, and Al₂O₃," *J. Opt. Soc. Am.* **65**, 742–744 (1975); detailed tables in DESY report SR-74/7 (1974).

No analytical Lorentz-Drude model exists for Mg (Rakić 1998 did not include it among its 11 metals). The Hagemann dataset covers an extraordinarily wide range: **0.0000248–24.8 µm** (50 meV to 500 keV), compiled from multiple measurement techniques via Kramers-Kronig analysis of reflectance combined with thin-film transmission measurements. The implementation here is **not** a full Hagemann table interpolation; it is a coarse piecewise approximation intended to capture the broad metallic behavior over the visible and infrared range.

The refractiveindex.info Mg entry lists three datasets: Hagemann et al. 1974 (the primary comprehensive set), Palm et al. 2018 (0.250–1.68 µm, more recent visible/near-UV data), and Vidal-Dasilva et al. 2010 (EUV/soft X-ray). Magnesium is **HCP and therefore optically anisotropic** (uniaxial); broad datasets generally represent a polycrystalline or effective average. Mg is highly reactive, forming a native MgO layer rapidly. For publication-grade Mg calculations, replacing the current coarse approximation with direct tabulated interpolation from Hagemann/Palm data would be preferable.

---

## 4. MgCl₂ uses a confirmed handbook value but lacks dispersion data entirely

**Value confirmed:** The refractive index **n = 1.675** corresponds to the ordinary-ray value (nω) of anhydrous MgCl₂ at the sodium D line (589 nm), consistent with the CRC Handbook and standard crystallographic databases. The extraordinary ray has nε = 1.59, giving a birefringence of δ = 0.085 that the code ignores.

Anhydrous MgCl₂ crystallizes in the trigonal system (space group R3̄m, CdCl₂ structure type). It is **extremely hygroscopic** – deliquescent to the point of dissolving in absorbed atmospheric moisture – which explains why no Sellmeier coefficients, detailed spectral dispersion data, or entries in Palik's *Handbook of Optical Constants* exist for this material. The refractiveindex.info database has no MgCl₂ entry.

The code's treatment (constant n = 1.675, k = 0 for λ < 2 µm, then piecewise IR absorption) is physically reasonable. The transparency window for a chloride ionic crystal likely extends from ~0.2 µm to somewhere in the 15–25 µm range by analogy with NaCl (0.21–26 µm) and KCl (0.21–30 µm). The onset of lattice absorption near 2 µm in the code seems conservative – real phonon absorption likely begins at longer wavelengths. However, without published single-crystal absorption spectra, the IR piecewise model **cannot be verified** against any known literature.

**Important caveat:** The hexahydrate MgCl₂·6H₂O (mineral: bischofite) has a significantly lower refractive index of ~1.569. The code value of 1.675 definitively refers to the anhydrous form.

---

## 5. MgF₂ coefficients from Dodge (1984) are exact, but C-value convention matters

**Source confirmed:** M. J. Dodge, "Refractive properties of magnesium fluoride," *Appl. Opt.* **23**, 1980–1985 (1984).

All six Sellmeier coefficients for the **ordinary ray** match the Dodge (1984) values exactly as listed on refractiveindex.info:

| Term | B | C (µm) |
|:----:|:-:|:------:|
| 1 | 0.48755108 | 0.04338408 |
| 2 | 0.39875031 | 0.09461442 |
| 3 | 2.3120353 | 23.793604 |

**Critical implementation detail:** The refractiveindex.info formula is n² − 1 = Σ Bᵢλ²/(λ² − Cᵢ²), where the C values are **resonance wavelengths that must be squared** in the denominator. They are not λ² values. The resonance wavelengths correspond to λ₁ = 43.4 nm (deep UV), λ₂ = 94.6 nm (VUV), and λ₃ = 23.8 µm (IR phonon). If the code treats C as already-squared (i.e., uses λ² − C instead of λ² − C²), the result at 589 nm would be n ≈ 1.439 instead of the correct **n = 1.378** – a significant error.

**Valid range: 0.2–7.0 µm** at 19°C, ordinary ray only. MgF₂ is uniaxial positive (tetragonal, rutile structure). The extraordinary ray requires separate Dodge-e coefficients. Alternative datasets include Li (1980) covering 0.14–7.5 µm and Zheng et al. (2023) providing temperature-dependent coefficients from 21–368°C.

---

## 6. Al₄C₃ has no published optical data – coefficients are unverified

**No literature source identified.** This is the only material in the set with coefficients that cannot be traced to any published reference. Extensive searching of refractiveindex.info, materials databases, and optical handbooks returned no Sellmeier equation, tabulated n,k data, or even isolated refractive index measurements for Al₄C₃.

The code's Sellmeier parameters (B₁ = 2.8, C₁ = 0.02; B₂ = 0.5, C₂ = 0.1; B₃ = 1.2, C₃ = 15.0) yield **n ≈ 2.08 at 589 nm** when C values are treated as resonance wavelengths. The previously reported n_D ≈ 2.7 from secondary sources could not be located in any accessible reference during this investigation. The discrepancy between 2.08 and 2.7 is ~30%, which is substantial for Mie calculations where scattering efficiency depends on |m − 1|².

What is known about Al₄C₃ optically comes primarily from Kioseoglou et al. (2019, *Phys. Status Solidi B* 256, 1900037), who measured an optical absorption bandgap of **~2.3 eV** (~539 nm) on single crystals. The material crystallizes in the trigonal system (R3̄m) and is therefore potentially birefringent, though the code uses a single isotropic formula. Al₄C₃ reacts violently with water (producing methane), making optical sample preparation extraordinarily difficult and explaining the absence of published data.

The resonance wavelengths in the code (20 nm, 100 nm, 15 µm) are physically plausible for two UV electronic transitions and one IR phonon, and the long-wavelength limit n∞ ≈ 2.35 is reasonable for a high-index ionic ceramic. However, the **reliability of these coefficients is low**, and users should treat results for Al₄C₃ particles with significant uncertainty.

---

## 7. MgO dispersion from Stephens & Malitson (1952) matches exactly; phonon model tracks Hofmeister (2003)

**Source confirmed:** R. E. Stephens and I. H. Malitson, "Index of refraction of magnesium oxide," *J. Res. Natl. Bur. Stand.* **49**, 249–252 (1952).

The transparent-region formula is a modified Ketteler-Helmholtz dispersion equation (not a pure Sellmeier), and all five coefficients match exactly:

$$n^2 = 2.956362 + \frac{0.02195770}{\lambda^2 - 0.01428322} - 0.01062387\lambda^2 - 2.04968 \times 10^{-5}\lambda^4$$

**Valid range: 0.36–5.4 µm** at 23.3°C. MgO is cubic (rocksalt structure, Fm3̄m) and therefore optically isotropic – no birefringence issues. The thermal coefficient dn/dT ranges from ~13.6 × 10⁻⁶/°C at 768 nm to ~19.0 × 10⁻⁶/°C at 405 nm.

For the **phonon region (λ > 5.4 µm)**, the code uses a four-oscillator classical dispersion model with ε∞ = 3.014 and oscillator centers at 384, 405, 429, and 590 cm⁻¹. These are consistent with Hofmeister, Keppel & Speck (2003, *MNRAS* 345, 16–38), who fit MgO reflectivity with a threefold TO multiplet near 400 cm⁻¹ plus multiphonon absorption near 590 cm⁻¹. The main TO mode at **405 cm⁻¹** is confirmed by THz spectroscopy (Han et al. 2008, 12.03 THz ≈ 401 cm⁻¹), the 384 and 429 cm⁻¹ components represent asymmetric broadening of the fundamental, and the **590 cm⁻¹** feature arises from two-phonon combination processes. The code's ε∞ = 3.014 is slightly lower than the most commonly cited value of **3.0176** (Fontanella et al. 1974), a minor discrepancy (~0.1%) that has negligible impact on calculated optical constants.

---

## 8. Al₂O₃ sapphire coefficients trace to Malitson & Dodge (1972), not Malitson (1962)

**Source confirmed – but the attribution needs correction.** The code coefficients match exactly with the **Malitson & Dodge (1972)** Sellmeier formula, later published in M. J. Dodge, "Refractive index," in *Handbook of Laser Science and Technology*, Vol. IV, CRC Press (1986). This is distinct from the original Malitson (1962) paper, which provided different ordinary-ray coefficients only.

All 12 coefficients (6 ordinary, 6 extraordinary) match the refractiveindex.info entries to every significant digit:

| Ray | B₁ | C₁ (µm) | B₂ | C₂ (µm) | B₃ | C₃ (µm) |
|:---:|:--:|:-------:|:--:|:-------:|:--:|:-------:|
| o | 1.4313493 | 0.0726631 | 0.65054713 | 0.1193242 | 5.3414021 | 18.028251 |
| e | 1.5039759 | 0.0740288 | 0.55069141 | 0.1216529 | 6.5927379 | 20.072248 |

As with MgF₂, the **C values are resonance wavelengths, not wavelength-squared**. The formula form is n² − 1 = Σ Bᵢλ²/(λ² − Cᵢ²). The resonance wavelengths correspond to two UV electronic transitions (~73 nm and ~119 nm) and one IR lattice resonance (~18–20 µm).

**Valid range: 0.20–5.0 µm** at 20°C. Sapphire (α-Al₂O₃) is uniaxial negative with birefringence Δn ≈ −0.008 at visible wavelengths. These coefficients apply to bulk single-crystal synthetic sapphire; thin-film or amorphous Al₂O₃ has significantly lower n (~1.6–1.7). Beyond ~5 µm, the Sellmeier formula loses accuracy as multiphonon absorption grows toward the reststrahlen band (12–25 µm). For extended IR coverage, tabulated data from Querry (1985, 0.21–55.6 µm) is needed.

---

## 9. Brass Cu70/Zn30 uses tabulated Querry (1985) data with linear interpolation

**Source confirmed:** M. R. Querry, "Optical constants," Contractor Report CRDC-CR-85034, U.S. Army Chemical Research, Development and Engineering Center (1985). Data available at [refractiveindex.info](https://refractiveindex.info/?shelf=other&book=Cu-Zn&page=Querry-Cu70Zn30).

**Mathematical model – tabulated n, k with linear interpolation (np.interp):**

The code stores 28 data points from the Querry dataset spanning **0.21–50 µm** and uses NumPy linear interpolation between them. This is the standard approach for metals where no compact analytical model (like Lorentz-Drude) has been fitted. The tabulated values are for the Cu70/Zn30 composition (cartridge brass, α-phase), which is the typical alloy used in brass pigments by manufacturers such as U.S. Bronze and Eckart.

**Selected values from the code's table (matching Querry 1985):**

| λ (µm) | n | k |
|-------:|------:|------:|
| 0.40 | 1.445 | 1.805 |
| 0.50 | 0.686 | 2.250 |
| 0.55 | 0.527 | 2.765 |
| 0.60 | 0.450 | 3.253 |
| 0.70 | 0.446 | 4.106 |
| 1.00 | 0.603 | 6.367 |
| 2.00 | 1.711 | 13.10 |
| 5.00 | 7.097 | 29.90 |
| 10.0 | 16.88 | 51.60 |
| 20.0 | 35.82 | 86.38 |
| 50.0 | 110.5 | 170.2 |

**Physical behavior:**
- **Visible range (0.4–0.7 µm):** n < 1, k ≈ 2–4. This is where the characteristic golden color of brass originates – strong absorption in the blue region and suppression of copper interband transitions near 0.5–0.6 µm.
- **Near/mid-IR (1–10 µm):** k rises steeply following the Drude free-electron model; n increases monotonically.
- **Far-IR (10–50 µm):** Both n and k are large and continue to grow – brass behaves as a good metallic reflector.

**Density:** 8530 kg/m³ (8.53 g/cm³) for Cu70/Zn30, the standard handbook value for cartridge brass (UNS C26000).

**Key caveats:** The Querry 1985 data is the only publicly available broadband measurement of brass optical constants. It was measured on bulk polished alloy samples. Real brass pigment particles (thin flakes produced by U.S. Bronze, Eckart, and others) may differ slightly from bulk values due to surface roughness, oxide films, and particle shape effects. The original Querry dataset has a gap at λ = 15 µm; the code's linear interpolation bridges this gap between the 10 and 20 µm data points. For wavelengths outside the 0.21–50 µm range, np.interp clamps to the boundary values. A second composition, Cu90/Zn10 (red brass/tombac), was also measured by Querry and is available on refractiveindex.info, but the code implements only Cu70/Zn30 as it is more representative of commercial brass pigments.

---

## 10. ZnCl₂ has no measured optical tables, but DFPT + IR/Raman constrain a plausible model

**No published n(λ),k(λ) data exist for solid ZnCl₂ in any optical database** – refractiveindex.info has no entry, Palik's *Handbook of Optical Constants* does not include it, and the only Querry et al. infrared paper on ZnCl₂ (*Appl. Opt.* **17**, 3587, 1978) treats *aqueous* ZnCl₂, dominated by water bands. The reason is the salt's deliquescence: anhydrous α-ZnCl₂ dissolves in atmospheric moisture within minutes, making single-crystal IR reflectance measurements extraordinarily difficult. Despite this, handbook refractive indices, DFPT dielectric constants, and vibrational spectroscopy provide enough constraints for a physically motivated reconstruction.

**Visible/UV – handbook value confirmed by DFPT.** The CRC Handbook (Weast 1988) lists n_o = 1.681, n_e = 1.713 at λ = 589 nm for α-ZnCl₂ (tetragonal, space group I-4̄2d, No. 122; cristobalite-like network of corner-sharing ZnCl₄ tetrahedra). Birefringence Δn = +0.032 makes the crystal **uniaxial positive**. These values are independently verified by density-functional perturbation theory in the Materials Project (entry **mp-22909**, Petousis et al., *Sci. Data* **4**, 160134, 2017), which gives a high-frequency dielectric tensor with ε∞_⊥ = 2.83 and ε∞_∥ = 2.96. The implied electronic-only indices are √ε∞_⊥ = **1.682** and √ε∞_∥ = **1.720** – matching the CRC values to 0.06% and 0.4% respectively. For polycrystalline particles in Mie calculations, the orientation average n_avg = √[(2 n_o² + n_e²)/3] = **1.692** is appropriate.

**Band gap and UV transparency.** The Materials Project PBE band gap is E_g ≈ 4.06 eV (~306 nm); accounting for the ~30–50% PBE underestimate, the true optical gap lies near 5–6 eV (~210–250 nm). Hu & Mackenzie (*Appl. Phys. Lett.* **33**, 57, 1978) report the visible-UV absorption edge of high-purity dehydrated ZnCl₂ glass at **210–220 nm**. At 0.4 µm (3.1 eV) the material is therefore well below the band edge and k can be set to zero to numerical precision.

**IR – four-oscillator Lorentz reconstruction.** No measured single-crystal reflectivity exists for α-ZnCl₂, but the Lyddane–Sachs–Teller framework can be constrained using independent inputs:

| Input | Value | Source |
|:------|:------|:-------|
| ε∞ (orientation average) | 2.87 | DFPT, mp-22909 |
| ε₀ (orientation average) | 5.24 | DFPT, mp-22909 |
| Δε_lattice = ε₀ − ε∞ | 2.37 | derived |
| TO mode positions (cm⁻¹) | 305, 230, ~110, ~60–80 | Angell & Wong 1970; Janz & James 1974; Kalampounias 2006 |

The dominant IR-active TO is the **antisymmetric Zn–Cl stretch at 305 cm⁻¹** (~33 µm reststrahlen), the symmetric stretch at 230 cm⁻¹ is primarily Raman-active with weaker IR oscillator strength, and the bend/network modes lie below 110 cm⁻¹. The code partitions the 2.37 lattice contribution among these mode groups to obtain the following approximate Lorentz parameters for ε(ω) = ε∞ + Σⱼ Δεⱼ ωⱼ²/(ωⱼ² − ω² − iωγⱼ):

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 305 | 32.8 | 1.60 | 25 | ν₃ Zn–Cl antisymmetric stretch (T₂) |
| 2 | 230 | 43.5 | 0.50 | 20 | ν₁ Zn–Cl symmetric stretch (A₁) |
| 3 | 105 | 95   | 0.20 | 30 | ν₂ Cl–Zn–Cl bend |
| 4 | 60  | 167  | 0.07 | 30 | network / acoustic-edge modes |

These reproduce ε₀ = ε∞ + Σ Δεⱼ = 2.87 + 2.37 = **5.24**, satisfying the static-dielectric sum by construction. **The parameters are not directly fitted to a measured reflectivity spectrum** – none exists – and should be treated as a physically constrained reconstruction, not a Palik-grade dataset. The strongest reststrahlen feature near 33 µm and the absence of significant absorption shortward of ~12 µm are qualitative predictions of this model.

**Multiphonon edge.** Hu & Mackenzie (1978) and Almeida (*J. Am. Ceram. Soc.* **72**, 537, 1989) place the multiphonon transmission edge of ZnCl₂ glass near **12–13 µm** (50% transmission, 1 mm sample) extending to ~20 µm at the 1 dB/cm level. For α-ZnCl₂ crystal the multiphonon onset is expected to lie just shortward of the strongest TO at 33 µm, with the mid-IR window 5–12 µm essentially transparent. The current implementation does not include a separate empirical multiphonon baseline; any such term would need to be added and documented explicitly if better absorption data become available.

**Density: 2907 kg/m³.** This is the X-ray crystallographic density of α-ZnCl₂ (Z = 4, lattice parameters of Brehler 1961; Yakel & Brynestad 1978) and is consistent with CRC, Sigma-Aldrich, and VWR product specifications. The DFPT-relaxed cell in mp-22909 gives 2760 kg/m³, ~5% lower as expected from PBEsol's typical volume overestimate. Bulk poured powder of commercial anhydrous ZnCl₂ has a much lower apparent density (~1500–1800 kg/m³) due to porosity and partial hydration; this is irrelevant for Mie scattering on dense crystallites but matters if modeling porous aggregates.

**Key caveats.**
- **Hygroscopicity and polymorphism.** ZnCl₂ has four anhydrous polymorphs (α, β, γ, δ) plus a hexahydrate ZnCl₂·6H₂O with significantly lower density (~1.77 g/cm³) and different refractive index. The values 1.681/1.713 and 2907 kg/m³ apply strictly to **anhydrous α-ZnCl₂**; reagent-grade material (~98%) is typically α-phase but absorbs water rapidly on exposure.
- **No measured IR reflectivity.** The Lorentz parameters above derive from DFPT (ε∞, ε₀) plus FTIR/Raman band positions on glassy and crystalline samples. They satisfy the LST sum rule but the partition of Δε_lattice among the four oscillators is a constrained reconstruction, not a fit. Users requiring publication-grade IR n,k for ZnCl₂ should commission a thin-film reflectance measurement under dry-N₂ purge.
- **Querry-1978 caution.** The aqueous-solution data of Querry, Cary & Waring (*Appl. Opt.* **17**, 3587, 1978) is sometimes cited as ZnCl₂ optical data; it is not – the spectrum is dominated by water (3, 6, 12–15 µm bands) modulated weakly by the Zn–Cl stretch near 290 cm⁻¹.

**Primary references:**
- C. A. Angell & J. Wong, *J. Chem. Phys.* **53**, 2053 (1970), DOI 10.1063/1.1674287.
- G. J. Janz & D. W. James, *Spectrochim. Acta A* **30**, 717 (1974), DOI 10.1016/0584-8539(74)80189-3.
- H. L. Yakel & J. Brynestad, *Inorg. Chem.* **17**, 3294 (1978), DOI 10.1021/ic50189a067.
- H. Hu & J. D. Mackenzie, *Appl. Phys. Lett.* **33**, 57 (1978), DOI 10.1063/1.90183.
- A. C. Almeida, *J. Am. Ceram. Soc.* **72**, 537 (1989), DOI 10.1111/j.1151-2916.1989.tb06032.x.
- A. Sen, M. N. Rao, R. Mittal & S. L. Chaplot, *J. Phys.: Condens. Matter* **17**, 6179 (2005), arXiv:cond-mat/0509329.
- A. G. Kalampounias, S. N. Yannopoulos & G. N. Papatheodorou, *J. Chem. Phys.* **124**, 014504 (2006), DOI 10.1063/1.2151888.
- G. Petousis et al., *Sci. Data* **4**, 160134 (2017); Materials Project entry mp-22909.
- R. C. Weast (ed.), *CRC Handbook of Chemistry and Physics*, 68th ed. (1988) – n_o, n_e, density.

---

## Conclusion

The implementation combines several well-traced literature models with approximate models where complete optical constants are unavailable. Aluminum, MgF₂, MgO (transparent region), Al₂O₃, and brass Cu70/Zn30 are tied to published coefficients or tabulated values. Carbon soot follows the widely cited Chang & Charalampopoulos (1990) polynomial, though the fitted curve slightly underpredicts some visible-band tabulated values. Magnesium is currently represented by a coarse piecewise approximation rather than a full tabulated interpolation. MgCl₂ correctly uses a handbook n_D value but lacks published dispersion or absorption data. ZnCl₂ likewise has no measured optical tables; its visible refractive indices are independently supported by DFPT (Materials Project mp-22909), while its IR dispersion is a physically constrained Lorentz reconstruction based on DFPT ε∞/ε₀ and FTIR/Raman phonon-band positions.

The least reliable entry is **Al₄C₃**, whose Sellmeier coefficients have no identifiable source and yield n ≈ 2.08 at 589 nm – potentially 30% below secondary-source estimates. This should be treated as a provisional placeholder rather than a validated optical-constants model. A second practical finding cuts across multiple materials: the Sellmeier C parameters for MgF₂, Al₂O₃, and Al₄C₃ are resonance wavelengths in µm that must be squared in the denominator (λ² − C²). Misinterpreting them as pre-squared values introduces systematic errors of several percent. Finally, the valid wavelength ranges cluster around **0.2–5 µm for the oxide and fluoride Sellmeier models**, with separate phonon or tabulated models required for mid-IR to far-IR Mie calculations beyond ~5 µm.
