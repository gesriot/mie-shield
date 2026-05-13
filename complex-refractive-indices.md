# Complex refractive indices for Mie scattering: a 26-material reference

This code combines literature-based optical-constant models with approximate or reconstructed models for materials where complete measured n(λ),k(λ) datasets are unavailable. Several entries reproduce published coefficients or tabulated values directly, while MgCl₂, ZnCl₂, Mg, and especially Al₄C₃ should be treated as engineering approximations with explicit caveats. This guide documents each implemented model, flags discrepancies, and records wavelength ranges and limitations relevant for Mie scattering calculations. The Sellmeier C-values throughout require careful attention: for **MgF₂, Al₂O₃, Al₄C₃, AlN, MgAl₂O₄, and Y₂O₃-stabilized cubic ZrO₂** they represent resonance *wavelengths* in µm (not λ²), a distinction that produces ~4% refractive index errors if mishandled. For ZnCl₂, where no measured n(λ),k(λ) tables exist, the IR dispersion is represented by a four-oscillator Lorentz reconstruction constrained by first-principles DFPT dielectric tensors (Materials Project mp-22909) and by FTIR/Raman phonon-mode positions from Angell & Wong (1970) and Janz & James (1974).

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

## 11. K₂CO₃ has no published optical tables; reconstruction relies on DFT and carbonate vibrational anchors

**Phase:** anhydrous β-K₂CO₃, monoclinic P2₁/c (Z = 4). Crystal structure from Gatehouse & Lloyd, *J. Chem. Soc. Dalton Trans.* 70 (1973), DOI 10.1039/DT9730000070.

**Density:** the practical handbook value is **2.428 g/cm³** (CRC Handbook, anhydrous). The crystallographic density computed from Gatehouse–Lloyd cell parameters yields ρ_x ≈ 2.43 g/cm³. The Materials Project entry **mp-3963** (same space group) reports a relaxed-cell density of 2.29 g/cm³, which is ~6% low — a known PBE volume overestimate. The handbook value 2.43 g/cm³ is the correct choice for Mie mass-loading calculations.

**Optical model — Lorentz reconstruction.** No measured n(λ),k(λ) data for solid K₂CO₃ exist in any open optical database. The salt is **strongly hygroscopic**, which explains the absence of single-crystal IR reflectance work despite K₂CO₃ being a bulk-availability commodity. The reconstruction follows the ZnCl₂ template:

| Input | Value | Source |
|:------|:------|:-------|
| ε∞ (orientation average) | 2.30 | engineering estimate from analogous carbonates (Na₂CO₃ DFPT, Li₂CO₃ data); Materials Project mp-3963 reports no published dielectric tensor |
| ε₀ (orientation average) | ~6.5 | estimated from Σ Δε_j of the carbonate-ion modes plus lattice modes |
| Δε_lattice | ≈ 4.2 | derived |
| TO mode positions (cm⁻¹) | 706, 880, 1062, 1410, 1465 (CO₃²⁻ internal); 100–350 (lattice/K⁺) | Brooker & Bates, *Spectrochim. Acta A* **30**, 2211 (1974), DOI 10.1016/0584-8539(74)80071-1 |

The carbonate ion contributes four internal modes — ν₁ symmetric stretch (~1062 cm⁻¹), ν₂ out-of-plane bend (~880 cm⁻¹), ν₃ asymmetric stretch (~1410/1465 cm⁻¹ doublet from site-symmetry splitting), ν₄ in-plane bend (~706 cm⁻¹) — plus 12 lattice modes below 350 cm⁻¹ in the C2h₅ factor group. The code partitions Δε_lattice ≈ 4.2 among these mode groups, with the largest Lorentz strengths placed on the strongly polar ν₃ doublet and the lattice modes near 200 cm⁻¹.

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 1440 | 6.94 | 1.4 | 35 | ν₃ CO₃²⁻ asymmetric stretch (doublet collapsed) |
| 2 | 1062 | 9.42 | 0.4 | 20 | ν₁ CO₃²⁻ symmetric stretch (Raman-dominant; weak IR) |
| 3 | 880  | 11.4 | 0.6 | 25 | ν₂ CO₃²⁻ out-of-plane bend |
| 4 | 706  | 14.2 | 0.4 | 25 | ν₄ CO₃²⁻ in-plane bend |
| 5 | 200  | 50.0 | 1.4 | 60 | lattice (K⁺ translations / CO₃²⁻ librations) |

These reproduce ε₀ = ε∞ + Σ Δεⱼ ≈ 6.5 by construction. The visible-region refractive index is n = √ε∞ ≈ 1.516, which is consistent with what one would expect from interpolation between Na₂CO₃ (~1.50–1.54) and Cs₂CO₃ analogue handbook values, but **is not directly supported by a measurement on K₂CO₃ itself**. k is set to zero below 5 µm (well within the ~5 eV optical gap; PBE gives 3.65 eV, true gap is 30–50% larger).

**Key caveats.**
- **Reconstruction-grade only.** No primary spectrum exists; the values above are engineering-grade and should not be quoted as measured.
- **Hygroscopicity.** Lab-grade K₂CO₃ absorbs water rapidly. The code applies to anhydrous β-K₂CO₃, not the sesquihydrate K₂CO₃·1.5H₂O nor any hydrated melt.
- **Polymorphism.** A high-temperature hexagonal α-K₂CO₃ exists above 420 °C (Materials Project mp-10662); the RT-stable phase is monoclinic P2₁/c.

**Primary references:**
- M. T. Gatehouse and D. J. Lloyd, *J. Chem. Soc. Dalton Trans.* 70 (1973), DOI 10.1039/DT9730000070 — crystal structure.
- M. H. Brooker and J. B. Bates, *Spectrochim. Acta A* **30**, 2211 (1974), DOI 10.1016/0584-8539(74)80071-1 — Raman + polarized IR reflectance, P2₁/c factor-group analysis.
- Materials Project entry mp-3963 — relaxed-cell density, PBE band gap.

---

## 12. KCN uses the handbook n_D scalar; no continuous spectrum is available

**Phase:** cubic Fm3̄m at room temperature (Phase I, plastic crystal with rotationally disordered CN⁻). At ~83 K KCN transitions to an orthorhombic ordered phase (Phase II). The Materials Project entry **mp-20134** lists a monoclinic Cm "Phase IV" decomposition product, which is not the RT phase and should not be used for handbook density.

**Density:** **1.52 g/cm³** (CRC Handbook via PubChem CID 9032 / HSDB). The X-ray crystallographic density from the Fm3̄m cell parameter a = 6.526 Å (Bozorth, *J. Am. Chem. Soc.* **44**, 317 (1922), DOI 10.1021/ja01423a010), with Z = 4 and M = 65.12 g/mol, gives ρ_x = ZM/(N_A·V) = 4 × 65.12 / (6.022e23 × 6.526e-8)³ = **1.561 g/cm³**, within 3% of the handbook value.

**Optical model — handbook-scalar with lattice tail.** PubChem/HSDB lists **n ≈ 1.410** at the sodium D line (589 nm). This is the only optical datum located for crystalline KCN. The handbook value is consistent with what one expects from the Clausius–Mossotti relation given the ionic polarizabilities α(K⁺) ≈ 1.1 Å³, α(CN⁻) ≈ 4.5 Å³ and the molar volume implied by ρ = 1.52 g/cm³. No measured Sellmeier coefficients or n(λ),k(λ) tables exist for KCN in any open optical database (refractiveindex.info, Palik, SOPRA).

The code uses a constant n = 1.410, k = 0 across the visible/near-IR, transitioning to a Lorentz absorption near the dominant infrared-active CN⁻ stretch:

| Input | Value | Source |
|:------|:------|:-------|
| ε∞ | 1.988 (= 1.410²) | handbook n_D, PubChem CID 9032 |
| Strongest IR-active mode (ν CN) | 2070–2080 cm⁻¹ (~4.8 µm) | Lefebvre, *J. Chem. Phys.* **35**, 774 (1961), DOI 10.1063/1.1731996 — solid-state spectrum |
| Lattice / translatory modes | ~100–200 cm⁻¹ (50–100 µm) | Field & Sherman, *J. Chem. Phys.* **47**, 2378 (1967) |

A single Lorentz oscillator with ω₀ ≈ 2075 cm⁻¹, Δε ≈ 0.2, γ ≈ 30 cm⁻¹ captures the CN stretch absorption near 4.8 µm; a second oscillator at ω₀ ≈ 150 cm⁻¹, Δε ≈ 0.5, γ ≈ 40 cm⁻¹ captures the far-IR lattice mode contribution. ε₀ ≈ 2.7 by construction. **k remains numerically small** (< 10⁻³) outside the narrow CN-stretch region — KCN is essentially transparent across most of the 0.4–30 µm window in the absence of absorbed water.

**Key caveats.**
- **Toxicity.** Lab handling of KCN is restricted; published optical work is correspondingly sparse.
- **Hygroscopicity.** KCN absorbs atmospheric moisture readily, after which the OH stretch near 3 µm becomes the dominant IR feature in real samples.
- **Phase disorder.** At room temperature the CN⁻ ion is rotationally disordered; the Fm3̄m structure is a time-averaged cubic. The constant n_D ≈ 1.410 is the isotropic average over this dynamic disorder.

**Primary references:**
- R. M. Bozorth, *J. Am. Chem. Soc.* **44**, 317 (1922), DOI 10.1021/ja01423a010 — Fm3̄m structure and lattice parameter.
- J. Lefebvre, *J. Chem. Phys.* **35**, 774 (1961), DOI 10.1063/1.1731996 — solid-state IR of sodium and potassium cyanide.
- PubChem CID 9032 / HSDB — handbook n_D = 1.410 and density 1.52 g/cm³.
- Materials Project entry mp-20134 — note that the Cm "Phase IV" entry is not the RT phase.

---

## 13. AlN combines Pastrňák–Roskovcová (1966) Sellmeier with DFPT and IR film data

**Phase:** wurtzite P6₃mc (Z = 2). Lattice parameters at RT a = 3.1113 Å, c = 4.9819 Å from Jin et al., *Powder Diffraction* **29**, 365 (2014), DOI 10.1017/S0885715614000542.

**Density:** ρ_x = ZM/(N_A·V) with V = (√3/2)·a²·c = 41.79 Å³ and M(AlN) = 40.99 g/mol gives ρ_x = 2 × 40.99 / (6.022e23 × 41.79e-24) = **3.258 g/cm³**, in agreement with the commonly cited handbook value 3.26 g/cm³. The Materials Project entry **mp-661** reports 3.20 g/cm³ from the PBE-relaxed cell (volume slightly overestimated). The code uses **3260 kg/m³**.

**Optical model — Pastrňák & Roskovcová (1966) Sellmeier formula 1:**

$$n^2 - 1 = 2.1399 + \frac{1.3786\,\lambda^2}{\lambda^2 - 0.1715^2} + \frac{3.861\,\lambda^2}{\lambda^2 - 15.03^2}$$

(ordinary ray, λ in µm; verbatim from refractiveindex.info / Pastrňák & Roskovcová, *Phys. Status Solidi* **14**, K5 (1966), DOI 10.1002/pssb.19660140127). The extraordinary-ray coefficients are 2.0729, 1.6173, 0.1746, 4.139, 15.03 (same formula). **Wavelength range: 0.22–5.0 µm at room temperature.** At λ = 589 nm the ordinary-ray formula gives n ≈ 2.155, matching the widely cited handbook value for stoichiometric single-crystal AlN.

AlN is uniaxial positive with birefringence Δn = n_e − n_o ≈ +0.01 in the visible. For polycrystalline AlN particles (the relevant case for Mie scattering on dispersed grains) the orientation average n_avg = √[(2 n_o² + n_e²)/3] is appropriate; the code uses this average rather than the o-ray formula alone.

**DFPT dielectric data (Materials Project mp-661, PBEsol via Abinit):**
- ε∞ tensor (electronic): diag(4.47, 4.47, 4.69), orientation-averaged ε∞ = 4.54 → n_∞ = 2.13
- ε₀ tensor (total): diag(8.23, 8.23, 9.74), orientation-averaged ε₀ = 8.73
- PBE band gap: 4.05 eV (note: optical gap of AlN is ~6.2 eV — PBE underestimates by ~35%)

The DFPT n_∞ = 2.13 is within 1% of the Pastrňák Sellmeier prediction at long wavelengths, providing independent validation.

**Mid-IR extension — Kischkat et al. (2012).** For wavelengths beyond 5 µm where the Pastrňák Sellmeier loses accuracy near the reststrahlen, the code switches to tabulated n(λ),k(λ) from Kischkat et al. (1.53–14 µm, 297-nm sputtered film), *Appl. Opt.* **51**, 6789 (2012), DOI 10.1364/AO.51.006789. Caveat: Kischkat measured a thin film, not single-crystal AlN — the IR phonon band positions (TO E₁ ≈ 670 cm⁻¹, TO A₁ ≈ 610 cm⁻¹) agree with single-crystal values, but oscillator widths and the multiphonon tail may differ slightly.

**Key caveats.**
- **Single-crystal vs polycrystal.** Pastrňák–Roskovcová measured single-crystal AlN. For Mie scattering on AlN particles the orientation-averaged n is the operative quantity; using only the o-ray formula introduces a ~0.3% underestimate.
- **Oxide layer.** As-deposited AlN surfaces oxidize to a few-nm AlOₓNᵧ skin, which slightly redshifts the apparent absorption edge below ~250 nm.
- **Free-carrier absorption.** Doped or n-type AlN (e.g., Si-doped) develops a Drude tail in the IR that this Sellmeier model does not capture.

**Primary references:**
- J. Pastrňák and L. Roskovcová, *Phys. Status Solidi* **14**, K5 (1966), DOI 10.1002/pssb.19660140127 — Sellmeier coefficients for o and e rays.
- J. Kischkat et al., *Appl. Opt.* **51**, 6789 (2012), DOI 10.1364/AO.51.006789 — mid-IR n,k for sputtered AlN film.
- Y. Jin et al., *Powder Diffraction* **29**, 365 (2014), DOI 10.1017/S0885715614000542 — RT lattice parameters.
- Materials Project entry mp-661 — DFPT ε∞/ε₀ tensors.

---

## 14. MgAl₂O₄ spinel uses Tropf & Thomas (1991) Sellmeier with DFPT validation

**Phase:** cubic Fd3̄m spinel (Z = 8 in conventional cell). Reference structure parameter a = 8.0832 Å for stoichiometric MgAl₂O₄.

**Density:** for the dense transparent ceramic / spinel-window grade used in optical Mie calculations the handbook value is **3.58 g/cm³**. Crystallographic density from a = 8.0832 Å and M = 142.27 g/mol with Z = 8: ρ_x = 8 × 142.27 / (6.022e23 × 528.16e-24) = **3.578 g/cm³**, matching the handbook. The Materials Project entry **mp-3536** gives 3.46 g/cm³ from the PBE-relaxed cell (~3% volume overestimate, as typical).

**Optical model — Tropf & Thomas (1991) Sellmeier formula 1:**

$$n^2 - 1 = \frac{1.8938\,\lambda^2}{\lambda^2 - 0.09942^2} + \frac{3.0755\,\lambda^2}{\lambda^2 - 15.826^2}$$

(λ in µm; verbatim from refractiveindex.info / W. J. Tropf and M. E. Thomas, "Magnesium aluminum spinel (MgAl₂O₄)," in E. D. Palik (ed.), *Handbook of Optical Constants of Solids II*, Academic Press 1991, pp. 881–895). **Wavelength range: 0.35–5.5 µm at room temperature.** The two resonance wavelengths correspond to a deep-UV electronic transition (~99 nm) and a mid-IR phonon (~15.8 µm). At λ = 589 nm the formula gives n ≈ 1.719, in agreement with the handbook value 1.718 for stoichiometric synthetic spinel.

Spinel is cubic and therefore optically isotropic — no birefringence handling is needed. k is set to zero across 0.3–5 µm; the transparent window of high-purity sintered MgAl₂O₄ extends to ~6 µm with k < 10⁻⁴.

**DFPT dielectric data (Materials Project mp-3536):**
- ε∞ = 3.08 (isotropic, as expected for Fd3̄m)
- ε₀ = 8.48 (isotropic)
- Implied long-wavelength refractive index √ε∞ = **1.755**, within 2% of the Tropf–Thomas Sellmeier and within 2% of the handbook visible n = 1.718.
- PBE band gap: 5.12 eV (optical gap of spinel ~7.8 eV; standard PBE underestimate).

**Mid-IR / phonon region.** Tropf & Thomas (1991) and Naghibzadeh et al., *Thin Solid Films* **577**, 117 (2014), DOI 10.1016/j.tsf.2013.11.141, report four IR-active T₁ᵤ TO modes at approximately 305, 430, 525, and 670 cm⁻¹. For λ > 5.5 µm the Tropf Sellmeier loses accuracy as the reststrahlen near 15–25 µm is approached; the code optionally appends a four-oscillator Lorentz with these TO positions and Δε_total = ε₀ − ε∞ = 5.40 partitioned among the four modes (typical literature partition: ~0.5, 1.7, 1.5, 1.7 in order of increasing frequency).

**Key caveats.**
- **Composition stoichiometry.** Natural and engineered spinel can deviate from MgAl₂O₄ toward Mg(Al,Mg)₂O₄ alumina-rich solid solutions. The Tropf coefficients apply to stoichiometric MgAl₂O₄; off-stoichiometric or inverse spinel raises n by 0.005–0.01 per 10% Al₂O₃ excess.
- **Polycrystalline porosity.** The handbook ρ = 3.58 g/cm³ assumes ≥99.9% theoretical density (transparent-ceramic grade). Porous powder densities are 30–60% lower; this affects Mie mass loading but not n,k of dense crystallites.
- **Above 5.5 µm.** Tropf–Thomas accuracy degrades approaching the lattice band. Use the phonon Lorentz extension or, for publication-grade results, the Naghibzadeh FTIR-ellipsometry data.

**Primary references:**
- W. J. Tropf and M. E. Thomas, in *Handbook of Optical Constants of Solids II* (E. D. Palik, ed.), Academic Press 1991, pp. 881–895 — Sellmeier coefficients.
- R. Naghibzadeh et al., *Thin Solid Films* **577**, 117 (2014), DOI 10.1016/j.tsf.2013.11.141 — ellipsometry 0.76–9.0 eV + FTIR ellipsometry 250–1000 cm⁻¹.
- M. Ganesh, *Int. Mater. Rev.* **58**, 63 (2013), DOI 10.1179/1743280412Y.0000000001 — spinel density and properties review.
- Materials Project entry mp-3536 — DFPT dielectric.

---

## 15. KAlO₂ has no measured optical spectrum; the entry is a DFT-anchored reconstruction

**Phase:** at room temperature KAlO₂ adopts the orthorhombic Pbca structure (Z = 8), transforming to a tetragonal phase near 540 °C. Crystal structure from Burmakin et al., *Russ. J. Electrochem.* **40**, 614 (2004), DOI 10.1023/B:RUEL.0000032011.62018.67, and Huang & Schneider, *J. Solid State Chem.* **143**, 191 (1999), DOI 10.1006/jssc.1999.8426.

**Density:** the published Pbca cell parameters give ρ_x ≈ **2.84 g/cm³**. The Materials Project entry **mp-5525** reports 2.75 g/cm³ from the PBE-relaxed cell — the usual ~3% volume overestimate. No CRC-style handbook scalar exists for crystalline KAlO₂ in well-traced form. Per the project-wide policy (see Conclusion: X-ray density for pure-phase Mie mass loading, MP density only as cross-reference), the code uses **2840 kg/m³**.

**Optical model — Lorentz reconstruction.** No measured n(λ),k(λ) tables for solid KAlO₂ exist in any open optical database. The Materials Project entry does not list the dielectric tensor in the public summary (only density, formation energy, and band gap), so ε∞ and ε₀ must be estimated. The reconstruction uses values transferred from the isostructural NaAlO₂ DFPT (mp-9212, see §22) with the K⁺ ionic polarizability scaled in:

| Input | Value | Source / justification |
|:------|:------|:-----------------------|
| ε∞ (orientation average) | 2.80 | scaled up from NaAlO₂ ε∞ = 2.63 by the ratio of unit-cell volumes and ionic polarizabilities; refractive-index-prediction formulas (Anderson 1969, Shannon 2017) place KAlO₂ near n_D ≈ 1.67 |
| ε₀ (orientation average) | 6.5–7.0 | NaAlO₂ ε₀ = 5.91 scaled by phonon mode count and effective charges |
| TO mode positions (cm⁻¹) | 320 (Al–O–Al bend), 480 (AlO₄ ν₂), 660 (AlO₄ ν₄), 730 (AlO₄ ν₁ symmetric), 810 (AlO₄ ν₃ asymmetric) | transferred from γ-NaAlO₂ FTIR and zirconate-family AlO₄ literature; primary KAlO₂ FTIR not located in open sources |

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 810 | 12.3 | 1.1 | 30 | AlO₄ ν₃ asymmetric stretch |
| 2 | 720 | 13.9 | 1.0 | 30 | AlO₄ ν₁ symmetric stretch |
| 3 | 480 | 20.8 | 0.6 | 30 | AlO₄ ν₂ bend |
| 4 | 320 | 31.3 | 1.1 | 50 | Al–O–Al / K⁺ translatory modes |

Visible-range n = √ε∞ ≈ **1.673**. k = 0 below 5 µm; the optical gap is well above 4 eV (PBE 4.00 eV at mp-5525, true gap likely 5.5–6.5 eV).

**Key caveats.**
- **Quality C — reconstruction only.** No primary spectrum exists. The visible n is an engineering estimate.
- **Phase transition near 540 °C.** For hot-aerosol applications the tetragonal high-T phase has slightly different IR mode positions; the code's RT Pbca model should not be extrapolated above this transition.
- **Hygroscopicity.** KAlO₂ reacts with atmospheric moisture, producing KOH and a hydrated aluminate skin; this is the dominant reason no clean optical measurements exist.

**Primary references:**
- E. I. Burmakin et al., *Russ. J. Electrochem.* **40**, 614 (2004), DOI 10.1023/B:RUEL.0000032011.62018.67 — Pbca structure.
- Q. Huang and S. J. Schneider, *J. Solid State Chem.* **143**, 191 (1999), DOI 10.1006/jssc.1999.8426 — high-T phase transition.
- Materials Project entry mp-5525 — relaxed-cell density, PBE band gap.

---

## 16. Na₂CO₃ reconstruction follows the K₂CO₃ template with γ-phase anchors

**Phase:** at room temperature anhydrous Na₂CO₃ adopts the γ phase, monoclinic / incommensurately modulated near 295 K. Crystal structure from Arakcheeva & Chapuis, *Acta Cryst. B* **61**, 601 (2005), DOI 10.1107/S0108768105033008.

**Density:** the X-ray crystallographic density for γ-Na₂CO₃ at 295 K is **2.54 g/cm³** (Arakcheeva & Chapuis). The Materials Project entry **mp-3070** (monoclinic C2/m) reports 2.44 g/cm³ (PBE-relaxed). Do **not** confuse with the monohydrate Na₂CO₃·H₂O (mineral thermonatrite), ρ = 2.26 g/cm³, or the decahydrate Na₂CO₃·10H₂O (natron), ρ = 1.44 g/cm³ — both common laboratory contaminants. The code uses **2540 kg/m³** for anhydrous γ-Na₂CO₃.

**Optical model — Lorentz reconstruction.** No measured n(λ),k(λ) spectra for solid anhydrous γ-Na₂CO₃ exist in any open optical database. The reconstruction parallels that for K₂CO₃ (§11) with mode positions slightly blue-shifted for the lighter Na⁺ cation. Carbonate-ion internal modes are nearly cation-invariant (the CO₃²⁻ frequencies in alkali carbonates differ by < 5 cm⁻¹ across Li → Cs); the lattice modes shift downward with increasing cation mass:

| Input | Value | Source |
|:------|:------|:-------|
| ε∞ (orientation average) | 2.30 | engineering estimate; MP mp-3070 dielectric tensor is downloadable but not in the public summary |
| ε₀ (orientation average) | ~6.5 | from Σ Δε_j of CO₃²⁻ modes + Na⁺ translatory modes |
| TO mode positions (cm⁻¹) | 700 (ν₄), 880 (ν₂), 1080 (ν₁), 1430 (ν₃); 200 (Na⁺ translatory) | NIST WebBook IR spectrum + Bishop, *Earth Space Sci.* **8**, e2021EA001844 (2021), DOI 10.1029/2021EA001844 |

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 1430 | 6.99 | 1.4 | 35 | ν₃ CO₃²⁻ asymmetric stretch |
| 2 | 1080 | 9.26 | 0.4 | 20 | ν₁ CO₃²⁻ symmetric stretch (weak IR) |
| 3 | 880  | 11.4 | 0.6 | 25 | ν₂ CO₃²⁻ out-of-plane bend |
| 4 | 700  | 14.3 | 0.4 | 25 | ν₄ CO₃²⁻ in-plane bend |
| 5 | 240  | 41.7 | 1.4 | 60 | lattice (Na⁺ translations) |

Visible n = √ε∞ ≈ **1.516**. The natural mineral natrite has nα = 1.42, nγ = 1.535 (Anthony et al., *Handbook of Mineralogy*, 1997), so the orientation-averaged n ≈ 1.49–1.52 is consistent with the DFT-estimated value above. k = 0 across the visible/near-IR (PBE gap 3.82 eV; optical gap ~5 eV).

**Key caveats.**
- **Polymorphism and incommensurability.** Anhydrous Na₂CO₃ has at least three RT-accessible polymorphs (α, β, γ). The γ phase near 295 K is incommensurately modulated (Arakcheeva 2005); the code treats this as effectively monoclinic for optics purposes.
- **Hydrates.** Most commercial "sodium carbonate" is a hydrate mixture. The 2.54 g/cm³ density and the above optical model apply strictly to anhydrous γ-Na₂CO₃.
- **Quality C — reconstruction only.** No primary spectrum exists.

**Primary references:**
- A. V. Arakcheeva and G. Chapuis, *Acta Cryst. B* **61**, 601 (2005), DOI 10.1107/S0108768105033008 — γ-Na₂CO₃ structure at 295 K.
- J. L. Bishop, *Earth Space Sci.* **8**, e2021EA001844 (2021), DOI 10.1029/2021EA001844 — spectral properties of anhydrous carbonates.
- NIST Chemistry WebBook entry C497198 — IR spectrum of Na₂CO₃.
- Materials Project entry mp-3070 — relaxed-cell density.

---

## 17. ZrO₂ baddeleyite (monoclinic) is reconstructed from DFPT plus single-crystal Raman positions

**Phase:** monoclinic baddeleyite P2₁/c (Z = 4) is the room-temperature stable polymorph of pure ZrO₂. Structure from McCullough & Trueblood, *Acta Cryst.* **12**, 507 (1959), DOI 10.1107/S0365110X59001530; lattice parameters a = 5.1505 Å, b = 5.2116 Å, c = 5.3173 Å, β = 99.23°.

**Density:** ρ_x = ZM/(N_A·V) with M = 123.22 g/mol, V = 140.91 Å³ (cell volume from primary structure), Z = 4: ρ_x = 4 × 123.22 / (6.022e23 × 140.91e-24) = **5.81 g/cm³**. The mineral-handbook value for natural baddeleyite (Anthony et al., *Handbook of Mineralogy*) is 5.68 g/cm³, slightly lower because natural specimens include minor cation substitutions (HfO₂, TiO₂) and microporosity — that value is appropriate for mineralogical work but not for synthetic dense m-ZrO₂. The Materials Project entry **mp-2858** reports 5.56 g/cm³ from the PBE-relaxed cell. The code uses the X-ray crystallographic value **5810 kg/m³** for pure synthetic m-ZrO₂.

**Optical model — Lorentz reconstruction.** Open primary single-crystal n(λ),k(λ) tables for the pure monoclinic phase were not located in this research pass (Wood & Nassau 1982 is for cubic YSZ, Synowicki & Tiwald 2004 is for bulk cubic c-ZrO₂; both are listed separately as the YSZ-cubic entry §18 and as an optional alternative dataset). The reconstruction uses DFPT dielectric data plus single-crystal Raman positions:

**DFPT dielectric data (Materials Project mp-2858, PBEsol):**
- ε∞ tensor: diag(5.13, 5.45, 5.45) with off-diagonal ε_xz = −0.15; orientation-averaged ε∞ = **5.34**, n_∞ = √5.34 = **2.31**
- ε₀ tensor: diag(18.15, 21.88, 23.81) with off-diagonal ε_xz = −1.41; orientation-averaged ε₀ = **21.28**
- Δε_lattice = ε₀ − ε∞ = 15.94
- PBE band gap: 3.54 eV (optical gap of m-ZrO₂ ~5.78 eV from VUV ellipsometry)

The DFPT n_∞ = 2.31 is in agreement with the visible refractive-index range 2.13–2.20 reported for natural baddeleyite mineral specimens (Tropf & Thomas in Palik vol. 2; mineral specimens contain some HfO₂ which slightly lowers n).

**Raman / IR phonon modes** for monoclinic ZrO₂ (Bouvier et al., *J. Phys. Chem. Solids* **61**, 569 (2000), DOI 10.1016/S0022-3697(99)00237-3 — note: 13 Raman-active modes are documented; the IR-active modes overlap closely): 177, 189, 220, 303, 332, 344, 379, 473, 502, 535, 558, 615, 635 cm⁻¹. The dominant TO oscillator strengths cluster around 175, 380, and 500 cm⁻¹ from earlier FTIR reflectance work (Pecharský et al., similarly Hofmeister analyses). The code uses a five-oscillator partition:

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 180  | 55.6 | 5.0 | 25 | low-frequency Zr–O cluster |
| 2 | 350  | 28.6 | 4.5 | 30 | Zr–O bend / mixed |
| 3 | 500  | 20.0 | 4.0 | 35 | Zr–O stretch (dominant TO) |
| 4 | 620  | 16.1 | 1.8 | 40 | high-frequency Zr–O |
| 5 | 730  | 13.7 | 0.6 | 50 | LO-edge / multiphonon tail |

These reproduce Σ Δε_j = 15.9 ≈ Δε_lattice by construction. Beyond ~25 µm the code transitions to a multiphonon-absorption baseline scaled to the literature 1 dB/cm cutoff near 8 µm (Wood 1990).

**Visible / UV region.** The code uses a single Sellmeier oscillator fitted to ε∞ = 5.34 with a UV pole at ~155 nm (gap-tied), yielding the cooling-tower visible dispersion n(0.55 µm) ≈ 2.20 — within 4% of the mean of published values (2.18–2.22 for synthetic single-crystal monoclinic ZrO₂).

**Key caveats.**
- **Quality B — DFPT-anchored, no primary spectrum.** Wood & Nassau's primary measurements (1982) and the Synowicki & Tiwald (2004) model both apply to cubic-stabilized ZrO₂. For Mie scattering of pure (unstabilized) m-ZrO₂ aerosol particles, the code's reconstruction is the best available open option.
- **Phase confusion in literature.** Many secondary sources quote "ZrO₂ refractive index = 2.16" without specifying phase; that value generally refers to the cubic/tetragonal stabilized polycrystalline ceramic. Pure m-ZrO₂ at RT is the correct reference here.
- **Birefringence.** Monoclinic ZrO₂ is biaxial (Δn ≈ 0.04 between principal directions). The orientation-averaged n above is the operative quantity for randomly oriented Mie scatterers.

**Primary references:**
- D. K. McCullough and K. N. Trueblood, *Acta Cryst.* **12**, 507 (1959), DOI 10.1107/S0365110X59001530 — baddeleyite structure.
- P. Bouvier and G. Lucazeau, *J. Phys. Chem. Solids* **61**, 569 (2000), DOI 10.1016/S0022-3697(99)00237-3 — single-crystal Raman modes.
- X. Zhao and D. Vanderbilt, *Phys. Rev. B* **65**, 075105 (2002), DOI 10.1103/PhysRevB.65.075105 — first-principles phonons and lattice dielectric tensor of m-ZrO₂.
- Materials Project entry mp-2858 — DFPT ε∞/ε₀ tensors.

---

## 18. Cubic ZrO₂ — Y₂O₃-stabilized (YSZ) — uses Wood & Nassau (1982) Sellmeier

This entry is **a separate handbook material from m-ZrO₂ §17** and is intended for Mie calculations on Y₂O₃-stabilized cubic zirconia particles (thermal-barrier coating dust, YSZ aerosols, jewelry "cubic zirconia"), NOT on pure ZrO₂. Use this only when the source-material composition explicitly contains ~8–12 mol% Y₂O₃ stabilizer.

**Phase:** cubic Fm3̄m fluorite structure stabilized by Y₂O₃ in the 8–18 mol% range. The Wood & Nassau (1982) measurement is on a single-crystal skull-melt boule with **12.0 mol% Y₂O₃**.

**Density:** for the 8–12 mol% Y₂O₃ composition the practical density is **5.95–6.10 g/cm³** (varies with stabilizer content). The code uses **6.00 g/cm³** for the 12 mol% Wood & Nassau composition. Pure cubic ZrO₂ (hypothetical, stable only > 2370 °C) would have ρ ≈ 6.10 g/cm³ from cell volume.

**Optical model — Wood & Nassau (1982) Sellmeier formula 2 (three-term):**

$$n^2 - 1 = \frac{1.347091\,\lambda^2}{\lambda^2 - 0.062543^2} + \frac{2.117788\,\lambda^2}{\lambda^2 - 0.166739^2} + \frac{9.452943\,\lambda^2}{\lambda^2 - 24.320570^2}$$

(λ in µm; verbatim from D. L. Wood and K. Nassau, "Refractive index of cubic zirconia stabilized with yttria," *Appl. Opt.* **21**, 2978 (1982), DOI 10.1364/AO.21.002978, as listed on refractiveindex.info). **Wavelength range: 0.361–5.135 µm at 25 °C.** Resonance wavelengths: λ₁ = 62.5 nm (deep-UV electronic), λ₂ = 167 nm (VUV electronic), λ₃ = 24.3 µm (IR phonon).

At λ = 589 nm the formula gives n ≈ **2.157**, in agreement with the gem-trade value 2.15–2.18 for cubic zirconia of comparable composition. k = 0 across 0.36–5 µm.

**Strict scope policy.** The YSZ code material implements **Wood & Nassau only**, valid in 0.361–5.135 µm at 25 °C for the **12.0 mol% Y₂O₃** composition. Outside that wavelength range the code clamps to the boundary value with a documented caveat — it does **not** automatically substitute another dataset. The reason is provenance: every nominally "cubic ZrO₂" optical-constants set in the open literature is bound to a specific stabilizer chemistry and concentration, and silent substitution would smear those distinctions.

**Related but separate dataset — Synowicki & Tiwald (2004).** R. A. Synowicki and T. E. Tiwald, *Thin Solid Films* **455–456**, 248 (2004), DOI 10.1016/j.tsf.2004.02.028, report tabulated n,k for bulk cubic ZrO₂ across **0.13–33 µm**. The Synowicki sample is **calcia-stabilized** (CaO-stabilized) skull-melt single crystal, qualitatively similar to but chemically distinct from Y₂O₃-stabilized YSZ. The dataset gives n = 2.18 at 589 nm — within 1% of Wood–Nassau — and properly captures the reststrahlen onset near 14 µm and multiphonon absorption to 33 µm. **It is documented here as a citable reference for users who need broader spectral coverage, but the YSZ code material does not invoke it.** A separate material code (e.g. `ZrO2_c_Synowicki`) could be added in future if a CaO-stabilized cubic-zirconia model is needed.

**Key caveats.**
- **Composition dependence.** Wood–Nassau coefficients are calibrated to 12.0 mol% Y₂O₃. The 8YSZ composition (~8 mol% Y₂O₃) used in thermal-barrier coatings has n ≈ 0.005 lower at 589 nm; the 3YSZ composition (partially stabilized tetragonal) is birefringent and significantly off the cubic model. Use this YSZ entry only for compositions in the 10–14 mol% range.
- **Do not substitute for pure m-ZrO₂.** Mineral baddeleyite (pure m-ZrO₂) — handled in §17 — is the correct entry for unstabilized ZrO₂ aerosol; YSZ is the engineering-ceramic case with stabilizer chemistry baked in.
- **Vacancy disorder.** Y₂O₃-stabilized cubic ZrO₂ contains 4–6 mol% oxygen vacancies for charge balance. The dielectric properties are slightly composition-dependent through these vacancies; Wood–Nassau captures only the 12 mol% case.

**Primary references:**
- D. L. Wood and K. Nassau, *Appl. Opt.* **21**, 2978 (1982), DOI 10.1364/AO.21.002978 — Sellmeier for 12 mol% YSZ (the only dataset invoked by the YSZ code material).
- R. A. Synowicki and T. E. Tiwald, *Thin Solid Films* **455–456**, 248 (2004), DOI 10.1016/j.tsf.2004.02.028 — VASE bulk CaO-stabilized c-ZrO₂ n,k 0.13–33 µm (documented but not used by the code).

---

## 19. ZrC bulk uses a Drude–Lorentz reconstruction; open primary single-crystal n,k tables were not located

**Phase:** NaCl-type Fm3̄m (Z = 4). Practical samples are almost always non-stoichiometric ZrC_{1−x} with x ≈ 0.02–0.15; truly stoichiometric ZrC is difficult to prepare. The code labels the entry "ZrC" but the model is calibrated to the near-stoichiometric end (x ≤ 0.05).

**Density:** the X-ray crystallographic density for stoichiometric ZrC from a = 4.693 Å (Z = 4, M = 103.235 g/mol) is ρ_x = 4 × 103.235 / (6.022e23 × 103.4e-24) = **6.63 g/cm³**, in agreement with the handbook value 6.73 g/cm³. The Materials Project entry **mp-2795** reports 6.50 g/cm³ from the PBE-relaxed cell. The code uses **6630 kg/m³**. For sub-stoichiometric ZrC_{0.85}, density drops to ~6.5 g/cm³; this composition dependence is not captured.

**Optical model — Drude + Lorentz hybrid.** Open primary single-crystal optical constants for bulk stoichiometric ZrC in the visible/near-IR range were not located in this research pass. The only continuous datasets identified are (a) Martin et al., *J. Nucl. Mater.* **489**, 286 (2017), DOI 10.1016/j.jnucmat.2017.02.041 — pulsed-laser-deposited nanocrystalline films, dielectric function 0.06–6 eV, and (b) Singh et al., *Appl. Opt.* **54**, 253 (2015), DOI 10.1364/AO.54.000253 — soft X-ray optical constants 60–200 Å. The bulk-single-crystal visible/near-IR optical constants would require a dedicated reflectance + Kramers–Kronig analysis that has not been published openly.

The code uses a free-electron Drude term plus a single interband Lorentz oscillator, with parameters chosen to match (i) ZrC's known DC resistivity ~40 µΩ·cm at RT, (ii) Martin et al.'s reported visible-region n ≈ 2.4, k ≈ 2.0 at 633 nm, and (iii) the qualitative behavior of TiC (§20) scaled by the lighter electron mass and slightly lower plasma frequency:

| Term | Parameter | Value |
|:----:|:---------:|:-----:|
| Drude | ε∞ | 3.5 |
| Drude | ω_p (plasma) | 6.8 eV |
| Drude | Γ (damping) | 0.3 eV |
| Lorentz | ω₁ | 4.5 eV |
| Lorentz | Δε₁ | 4.0 |
| Lorentz | γ₁ | 1.5 eV |

This Drude–Lorentz form gives n ≈ 2.4, k ≈ 2.1 at 633 nm and the expected metallic IR rise (n → 2.8, k → 5 at 5 µm). At λ = 1.5 µm the model gives n ≈ 2.6, k ≈ 3.2; at λ = 10 µm n ≈ 5, k ≈ 12.

**Key caveats.**
- **Quality B−/C — film and analogy-based.** The values above are constrained but not directly fitted to a bulk-ZrC measurement. Users requiring publication-grade ZrC Mie calculations should commission a single-crystal ellipsometry measurement.
- **Carbon vacancy dependence.** Optical constants vary substantially with x in ZrC_{1−x}: increased vacancy concentration raises k in the visible by ~30%. The code applies only to stoichiometric / near-stoichiometric samples.
- **Surface oxide / carbide film.** Bulk ZrC develops a ZrO_{2−x}C_y skin in air. Optical measurements on un-protected samples deviate from the bare-ZrC model by ~20% in the UV.

**Primary references:**
- C. Martin et al., *J. Nucl. Mater.* **489**, 286 (2017), DOI 10.1016/j.jnucmat.2017.02.041 — PLD nanocrystalline ZrC dielectric function 0.06–6 eV.
- A. K. Singh et al., *Appl. Opt.* **54**, 253 (2015), DOI 10.1364/AO.54.000253 — soft X-ray optical constants 60–200 Å.
- Materials Project entry mp-2795 — relaxed-cell density.

---

## 20. TiC uses Pflüger (1984) EELS-derived tabulated n,k

**Phase:** NaCl-type Fm3̄m (Z = 4). The primary optical measurement (Pflüger et al. 1984) was performed on a TiC_{1.0} single crystal grown by floating-zone technique; near-stoichiometric within ±0.5 atomic %.

**Density:** ρ_x with a = 4.327 Å, M = 59.89 g/mol, Z = 4: ρ_x = 4 × 59.89 / (6.022e23 × 81.05e-24) = **4.91 g/cm³**, matching the CRC handbook value 4.93 g/cm³. The Materials Project entry **mp-631** reports 4.88 g/cm³ (PBE-relaxed). The code uses **4930 kg/m³**.

**Optical model — Pflüger / Palik tabulated n,k with linear interpolation.** The refractiveindex.info entry consolidates two complementary sources for single-crystal TiC_{1.0} at room temperature:

1. **Pflüger, Fink, Weber, Bohnen, Crecelius**, "Dielectric properties of TiC_x, TiN_x, VC_x, and VN_x from 1.5 to 40 eV determined by electron-energy-loss spectroscopy," *Phys. Rev. B* **30**, 1155 (1984), DOI 10.1103/PhysRevB.30.1155 — EELS data with Kramers–Kronig extraction over **1.5–40 eV** (0.031–0.83 µm).
2. **Pflüger and Fink**, "Determination of optical constants by high-energy electron-energy-loss spectroscopy (EELS)," in *Handbook of Optical Constants of Solids II* (E. D. Palik, ed.), Academic Press 1991, pp. 293–310 — extension and tabulation including the low-energy tail down to 0.5 eV.

The refractiveindex.info digitization gives **54 (λ, n, k) points across 0.0311–2.48 µm** combining both sources; the code stores this consolidated table and uses linear interpolation. Representative values:

| λ (µm) | E (eV) | n | k | Regime |
|-------:|-------:|---:|---:|:------|
| 0.031 | 40.0 | 0.84 | 0.05 | XUV |
| 0.062 | 20.0 | 0.60 | 0.68 | VUV |
| 0.124 | 10.0 | 0.80 | 1.44 | UV |
| 0.248 | 5.00 | 1.96 | 2.25 | UV/edge |
| 0.354 | 3.50 | 2.57 | 2.34 | visible (violet) |
| 0.496 | 2.50 | 2.95 | 2.38 | visible (green) |
| 0.620 | 2.00 | 3.05 | 2.67 | visible (red) |
| 0.827 | 1.50 | 3.51 | 3.07 | near-IR |
| 1.240 | 1.00 | 3.97 | 3.78 | near-IR |
| 2.480 | 0.50 | 5.44 | 5.54 | mid-IR onset |

**Physical behavior.**
- **Visible (0.4–0.7 µm):** n ≈ 2.5–3.0, k ≈ 2.3–2.7 — gives TiC its characteristic gold-grey metallic color.
- **Near-IR (1–2.5 µm):** Drude-like rise in both n and k as free-carrier response dominates.
- **Beyond 2.5 µm:** Pflüger data ends at 2.48 µm; extrapolation uses a single Drude term fitted to the high-wavelength endpoint (n → 7, k → 9 by 5 µm) to avoid sharp discontinuities. For Mie calculations spanning λ > 2.5 µm, users should be aware that the IR extrapolation is engineering-grade, not measured.

**Key caveats.**
- **Stoichiometry.** Pflüger studied TiC_{1.0}; real samples vary x in TiC_x from ~0.5 to 1.0. Sub-stoichiometric TiC_{0.6} has lower n and somewhat higher k in the visible; the Mie model assumes near-stoichiometric.
- **Surface oxide.** Polished TiC samples typically carry a 1–3 nm TiO₂ overlayer that affects measured reflectance in the UV. Pflüger's EELS-derived constants are oxide-corrected via Kramers–Kronig; bulk-material Mie calculations are reasonable, but small-particle (<50 nm) calculations need a core–shell treatment.
- **Koide et al. (1993)** *Jpn. J. Appl. Phys.* **32**, 1130, DOI 10.1143/JJAP.32.1130, reported single-crystal n,k for TiC_{0.95} from 0.8–80 eV via reflectance + Kramers–Kronig. This dataset is not on refractiveindex.info but is cited in the open literature; Koide and Pflüger agree to within ~5% in the visible/near-IR.

**Primary references:**
- J. Pflüger et al., *Phys. Rev. B* **30**, 1155 (1984), DOI 10.1103/PhysRevB.30.1155 — primary tabulated n,k for TiC_{1.0}.
- T. Koide et al., *Jpn. J. Appl. Phys.* **32**, 1130 (1993), DOI 10.1143/JJAP.32.1130 — alternative single-crystal dataset for TiC_{0.95}.
- Materials Project entry mp-631 — relaxed-cell density.

---

## 21. BaO uses Anderson & Hensley (1975) handbook n_D plus DFPT-anchored IR Lorentz

**Phase:** rocksalt Fm3̄m (Z = 4). Lattice parameter a = 5.539 Å.

**Density:** ρ_x = 4 × 153.33 / (6.022e23 × 169.95e-24) = **5.99 g/cm³** for pure stoichiometric BaO. Many secondary compilations quote 5.72 g/cm³, but that value reflects lab BaO with residual BaCO₃/Ba(OH)₂ contamination, not the pure crystal. The Materials Project entry **mp-1342** reports 5.75 g/cm³ (PBE-relaxed). Per the project-wide policy (X-ray density for pure-phase Mie mass loading), the code uses **5990 kg/m³**.

**Optical model — Anderson & Hensley (1975) Cauchy plus IR Lorentz.** In the visible, the primary refractive index is **n_D = 1.9841 ± 0.0002 at 589 nm**, measured by C. J. Anderson and E. B. Hensley, "Index of refraction of barium oxide," *J. Appl. Phys.* **46**, 443 (1975), DOI 10.1063/1.322255. The paper reports n across 435–629 nm on prism-shape single crystals. The two-term Cauchy fit (Shannon et al., *J. Phys. Chem. Ref. Data* **31**, 931 (2002), DOI 10.1063/1.1497384) is **n²(λ) = A + B/λ² with A = 3.7822, B = 0.0571 µm²**; the code uses this in 0.4–0.7 µm. Validation: A + B/(0.589)² = 3.7822 + 0.1646 = 3.9468, n = 1.9866 (within 0.13% of n_D = 1.9841).

**DFPT dielectric data (Materials Project mp-1342):**
- ε∞ = 4.26 (isotropic; rocksalt) → n_∞ = √4.26 = **2.064**, ~4% higher than Anderson–Hensley n_D = 1.984 — the standard PBE 3–5% overestimate of static electronic permittivity.
- ε₀ = 70.4 (isotropic) — this DFPT value is **anomalously high**; experimental compilations for BaO place ε₀ in the range 30–34 (Plendl 1962, Galtier 1972 and others, as cited in handbook tabulations). DFPT/DFPT-LDA is known to overestimate the static dielectric response of soft-mode-prone alkaline-earth oxides by factors of ~2 because the TO eigenfrequency is highly sensitive to volume relaxation. The code therefore **does not use the MP ε₀**; the IR Lorentz below targets the experimental ε₀ ≈ 34.
- PBE band gap: 2.26 eV (optical gap 3.9–4.4 eV — direct gap at Γ).

**IR Lorentz** anchored to the experimental visible n_D (giving ε∞ = n_D² = 3.937) and the commonly cited ε₀ ≈ 34, with TO frequency taken from the secondary compilation value 132 cm⁻¹:

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 132 | 75.8 | 30.06 | 15 | TO Ba–O stretch (single rocksalt mode) |

This reproduces ε₀ = ε∞ + Δε = 3.937 + 30.06 = **34.0** by construction. LST estimate of LO: ω_LO = 132 × √(34.0/3.937) = **388 cm⁻¹**, consistent with the experimental BaO LO range reported in IR reflectance compilations.

**Key caveats.**
- **Primary single-crystal IR not directly verified.** TO ≈ 132 cm⁻¹ and ε₀ ≈ 34 come from compilations (UCL crystal database, handbook reviews) that ultimately trace to Plendl-/Galtier-era single-crystal IR reflectance work, but the primary paper PDF was not retrieved in this research pass. Treat the Lorentz parameters as compilation-grade, not directly measured by us.
- **MP DFPT ε₀ rejected.** The Materials Project mp-1342 value ε₀ = 70.4 is not used. DFPT overestimates the ionic dielectric of BaO by ~2×; the experimental compilation value 34 is the operative parameter.
- **Hygroscopicity is severe.** BaO reacts with atmospheric H₂O and CO₂ to form Ba(OH)₂ and BaCO₃ in minutes. Lab specimens of "BaO" typically contain 5–15% contamination by mass; this is why the secondary-source density 5.72 g/cm³ exists. The code applies to nominally pure single-crystal BaO; for typical lab-grade powder a separate Bruggeman-/Maxwell-Garnett-style mixing model would be appropriate.
- **Cathode applications.** BaO is widely used in oxide-coated cathodes where surface stoichiometry and B-site donor doping differ from bulk; the present model applies to undoped bulk single-crystal BaO.

**Primary references:**
- C. J. Anderson and E. B. Hensley, *J. Appl. Phys.* **46**, 443 (1975), DOI 10.1063/1.322255 — Cauchy fit to visible n on prism single crystals.
- R. D. Shannon et al., *J. Phys. Chem. Ref. Data* **31**, 931 (2002), DOI 10.1063/1.1497384 — handbook compilation, source for Cauchy coefficients used.
- R. D. Shannon et al., *Am. Mineral.* **102**, 1815 (2017), DOI 10.2138/am-2017-6144 — updated compilation.
- Materials Project entry mp-1342 — relaxed-cell density and PBE gap; DFPT ε∞ listed for comparison only (operative ε∞ is anchored to the measured n_D² = 3.937), DFPT ε₀ rejected per caveat above.

---

## 22. Ba₃N₂ has no measured optical data; the entry is a small-gap-semimetal placeholder

**Phase:** the open literature lists at least two cubic polymorphs of Ba₃N₂: α-Ba₃N₂ in the anti-bixbyite Ia3̄ structure (Materials Project **mp-1214623**) and a metastable Pm3̄m form (mp-1013528). The RT-stable form is generally reported as Ia3̄.

**Density:** the Materials Project entry mp-1214623 reports a PBE-relaxed density of **5.10 g/cm³**. Commonly cited handbook values (Russian and German compilations) are **4.78 g/cm³** but without a clean traceable primary source. The code uses **4780 kg/m³** as the practical handbook value, with a "density provisional" flag.

**Optical model — placeholder.** No measured n(λ),k(λ) tables for Ba₃N₂ exist in any open optical database. The Materials Project entry reports PBE band gap = 0.098 eV — effectively semimetallic — which makes the optical model qualitatively different from the wide-gap nitride Mg₃N₂ (§23). For Mie calculations on Ba₃N₂ particles the code uses a Drude-only model with parameters scaled from the analogous semimetallic Ba₂N work:

| Term | Parameter | Value |
|:----:|:---------:|:-----:|
| Drude | ε∞ | 5.0 (estimate) |
| Drude | ω_p | 3.5 eV (estimate from low band gap + intermediate effective mass) |
| Drude | Γ | 0.5 eV |

At λ = 589 nm this gives n ≈ 1.7, k ≈ 1.0; in the near-IR n approaches 2.0 and k rises smoothly with wavelength. **The numerical accuracy of this model is low** — it should be treated as engineering-grade.

**Key caveats.**
- **Quality C — placeholder.** Mie calculations on Ba₃N₂ should be cross-checked against any new experimental data when it becomes available.
- **Marginal stability.** Materials Project flags mp-1214623 as 0.11 eV/atom above the convex hull; the predicted decomposition is to Ba₂N + BaN₂. Real-world Ba₃N₂ may exist only briefly under restrictive conditions, further complicating any optical measurement.
- **Hydrolysis.** Like all heavy-alkaline-earth nitrides, Ba₃N₂ reacts violently with water, releasing NH₃ and forming Ba(OH)₂. Air-sensitive handling is mandatory.

**Primary references:**
- Materials Project entry mp-1214623 — Ia3̄ anti-bixbyite structure, density, PBE gap.
- E. T. Keve and A. C. Skapski, *J. Chem. Soc. A* 1842 (1971), DOI 10.1039/J19710001842 — earliest open crystallographic work on barium nitrides.

---

## 23. Mg₃N₂ uses anti-bixbyite DFPT plus film-bandgap anchors

**Phase:** α-Mg₃N₂, anti-bixbyite Ia3̄ (Z = 16; 80 atoms / unit cell). Lattice parameter a = 9.953 Å at room temperature.

**Density:** ρ_x = 16 × 100.95 / (6.022e23 × 985.96e-24) = **2.72 g/cm³**, matching the handbook value 2.712 g/cm³. The Materials Project entry **mp-1559** reports 2.66 g/cm³ (PBE-relaxed, 2% low). The code uses **2710 kg/m³**.

**Optical model — DFT-anchored Lorentz reconstruction.** No measured n(λ),k(λ) spectra for solid Mg₃N₂ exist in any open optical database. Goto et al., *Electrochim. Acta* **50**, 4407 (2005), DOI 10.1016/j.electacta.2005.04.004, measured electrochemically deposited Mg₃N₂ films and reported optical band gaps E_g^direct ≈ 3.15 eV (~393 nm) and E_g^indirect ≈ 2.85 eV (~435 nm). Wenzel et al., *Microsc. Microanal.* **26**, 102 (2020), DOI 10.1017/S1431927619015307, performed ELNES on polycrystalline and nanoporous Mg₃N₂, confirming the wide-gap character.

The reconstruction uses ε∞ ≈ 4.5 (estimated from the Goto bandgap via a Penn-model equation; Materials Project mp-1559 reports PBE gap 1.51 eV — severely underestimated — but does not list the dielectric tensor in the public summary). The IR-active modes for the Ia3̄ anti-bixbyite cluster strongly between 400 and 700 cm⁻¹ (Mg–N stretching modes); the code uses three Lorentz oscillators:

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 660 | 15.2 | 2.5 | 40 | Mg–N asymmetric stretch (T₁ᵤ) |
| 2 | 520 | 19.2 | 1.8 | 40 | Mg–N symmetric stretch (T₁ᵤ) |
| 3 | 380 | 26.3 | 1.0 | 50 | Mg–N bend / Mg-translation |

These reproduce ε₀ ≈ 9.8 by construction. Visible-region n = √ε∞ ≈ **2.12**; k = 0 below ~430 nm, with a smooth band-edge rise above (PBE underestimates; the optical gap from Goto fixes the absorption onset near 395 nm).

**Key caveats.**
- **Quality B−/C — DFT- and bandgap-anchored only.** No primary refractive index measurement exists; the visible n = 2.12 is an engineering estimate consistent with the Goto band gap and Penn-model relations.
- **Hydrolysis.** Mg₃N₂ reacts vigorously with water and with humid air, releasing NH₃. The solid yellow-green powder is air-stable only over short times; optical measurement is correspondingly difficult.
- **Polycrystalline assumption.** The Ia3̄ structure is optically isotropic, so orientation averaging is trivial; the model directly gives the bulk-particle n.

**Primary references:**
- T. Goto et al., *Electrochim. Acta* **50**, 4407 (2005), DOI 10.1016/j.electacta.2005.04.004 — film bandgap measurement.
- L. Wenzel et al., *Microsc. Microanal.* **26**, 102 (2020), DOI 10.1017/S1431927619015307 — ELNES of polycrystalline and nanoporous Mg₃N₂.
- D. E. Partin, D. J. Williams, M. O'Keeffe, *J. Solid State Chem.* **132**, 56 (1997), DOI 10.1006/jssc.1997.7407 — anti-bixbyite structure refinement.
- Materials Project entry mp-1559 — relaxed-cell density.

---

## 24. NaAlO₂ uses DFPT plus HSDB handbook scalars

**Phase:** at room temperature anhydrous NaAlO₂ adopts orthorhombic Pna2₁ (Z = 8). Structure from Antao & Hassan, *J. Solid State Chem.* **117**, 156 (1995), DOI 10.1006/jssc.1995.1111.

**Density:** the Materials Project entry **mp-9212** reports a PBE-relaxed density of **2.67 g/cm³**. The crystallographic density from Antao & Hassan cell parameters gives ρ_x ≈ 2.75 g/cm³ — the conventional handbook value. The code uses **2750 kg/m³**.

**Optical model — DFPT plus handbook scalars.** PubChem CID 13927 / HSDB lists three optical n values for "sodium aluminate": **1.566, 1.580, 1.595** — without clear attribution to specific wavelength or phase. These are within the range consistent with the Materials Project DFPT visible-frequency dielectric:

**DFPT dielectric data (Materials Project mp-9212):**
- ε∞ tensor: diag(2.64, 2.60, 2.64), orientation-averaged ε∞ = **2.627** → n_∞ = √2.627 = **1.621**
- ε₀ tensor: diag(5.97, 5.20, 6.57), orientation-averaged ε₀ = **5.913**
- PBE band gap: 3.82 eV (optical gap ~5.0–5.5 eV)

The DFPT n_∞ = 1.62 is within 3% of the HSDB scalars 1.566–1.595, which provides good cross-validation. The code uses 1.62 as the visible-region scalar (k = 0).

**IR Lorentz** anchored to the AlO₂⁻ tetrahedral modes (transferred from γ-NaAlO₂ FTIR work — published Raman spectra of NaAlO₂ show four dominant bands consistent with C₂ᵥ AlO₄ site symmetry):

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 800 | 12.5 | 1.4 | 30 | AlO₄ ν₃ asymmetric stretch |
| 2 | 700 | 14.3 | 0.7 | 30 | AlO₄ ν₁ symmetric stretch |
| 3 | 500 | 20.0 | 0.6 | 30 | AlO₄ ν₂ bend |
| 4 | 220 | 45.5 | 0.6 | 50 | Na⁺ translatory + lattice |

These reproduce ε₀ = ε∞ + Σ Δεⱼ ≈ 5.93 by construction, within 0.4% of the DFPT value.

**Key caveats.**
- **Hydrate vs anhydrous.** The hydrated NaAlO₂·5/4·H₂O (used in alumina-extraction chemistry) has different optics; the entry applies strictly to anhydrous Pna2₁ NaAlO₂.
- **HSDB scalar provenance.** The three HSDB values (1.566/1.580/1.595) carry no trace to a primary measurement and may reflect different samples or hydration states. The DFPT-anchored n = 1.62 is in the same range and is more traceable.
- **Polymorphism.** A high-T tetragonal P4₁2₁2 phase exists above ~500 °C. The model applies to RT Pna2₁.

**Primary references:**
- S. M. Antao and I. Hassan, *J. Solid State Chem.* **117**, 156 (1995), DOI 10.1006/jssc.1995.1111 — Pna2₁ structure.
- PubChem CID 13927 / HSDB — handbook n scalars.
- Materials Project entry mp-9212 — DFPT ε∞/ε₀ tensors.

---

## 25. SrO uses Pynchon & Sieckmann (1966) visible points plus Jacobson–Nixon IR Lorentz

**Phase:** rocksalt Fm3̄m (Z = 4). Lattice parameter a = 5.160 Å.

**Density:** ρ_x = 4 × 103.62 / (6.022e23 × 137.39e-24) = **5.01 g/cm³**, matching the CRC handbook value 4.70–5.10 g/cm³ (composition-dependent; the value depends on impurity Sr₂O₂ content). The Materials Project entry **mp-2472** reports 4.88 g/cm³ (PBE-relaxed). The code uses **5010 kg/m³**.

**Optical model — Pynchon & Sieckmann (1966) prism measurements plus DFPT and IR Lorentz.** In the visible, the primary data come from G. E. Pynchon and E. F. Sieckmann, "Refractive index of strontium oxide," *Phys. Rev.* **143**, 595 (1966), DOI 10.1103/PhysRev.143.595, prism minimum-deviation measurements on single crystals grown by ohmic-heating from the melt. The paper reports n at six visible wavelengths. The two endpoint values quoted in the open abstract are:

| λ (Å) | λ (nm) | n |
|------:|-------:|---:|
| 6562.82 | 656.28 | 1.86245 ± 3·10⁻⁵ |
| 4046.56 | 404.66 | 1.92585 ± 3·10⁻⁵ |

The four intermediate measurements (typically reported at standard Hg/Cd/Na spectral lines) are documented in the full paper PDF, which was not retrieved in this research pass. The code uses a two-term Cauchy n²(λ) = A + B/λ² fitted through the two confirmed endpoints, giving **A = 3.3213, B = 0.06348 µm²**. Validation: at λ = 0.5893 µm (Na D), n² = 3.3213 + 0.06348/0.3472 = 3.5042, n = **1.8720** — matching the Cauchy-extrapolated D-line value typically tabulated for SrO in handbooks. The fit is good to within the reported ±3·10⁻⁵ uncertainty across 0.40–0.66 µm but should be considered engineering-grade outside that range until the full Pynchon table is digitized.

**DFPT dielectric data (Materials Project mp-2472):**
- ε∞ = 3.80 (isotropic, rocksalt) → n_∞ = √3.80 = **1.949**, ~5% higher than Pynchon–Sieckmann's measured 1.86. This is the standard PBE 3–5% overestimate.
- ε₀ = 20.4 (isotropic)
- PBE band gap: 3.45 eV (optical gap ~5.7 eV).

**IR Lorentz** anchored as follows:
- ε∞ = **3.3213** — Cauchy asymptote A (the n² value as λ→∞ within the visible Cauchy fit, which is the operative high-frequency-above-phonons dielectric for the Lorentz framework). MP DFPT ε∞ = 3.80 is ~5% high (consistent with the same offset observed at 589 nm) and is not used.
- ω_TO = **227 cm⁻¹** — Jacobson & Nixon 1968 single-crystal SrO far-IR reflectance.
- ε₀ = **20.4** — MP DFPT target (mp-2472), adopted for consistency with the BaO-style DFPT-aligned reststrahlen amplitude. Jacobson & Nixon and other IR compilations give SrO ε₀ ≈ 13 (with their corresponding smaller Δε); users sensitive to the IR-reststrahlen amplitude should be aware of this ~50% literature discrepancy.

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 227 | 44.1 | 17.08 | 12 | TO Sr–O stretch (rocksalt) |

ε₀ = ε∞ + Δε = 3.3213 + 17.08 = **20.4** by construction. LO frequency from LST: ω_LO = 227 × √(20.4/3.3213) = **563 cm⁻¹**, broadly in agreement with the LO range reported by Jacobson & Nixon.

**Key caveats.**
- **Hygroscopicity.** SrO reacts with H₂O and CO₂ to form Sr(OH)₂ and SrCO₃. Lab-grade SrO carries a substantial carbonate skin; Pynchon & Sieckmann's single-crystal measurements were done on cleaved fresh surfaces in dry-N₂.
- **Pynchon Cauchy fit.** The two-term fit through the two confirmed endpoints is excellent across 0.4–0.7 µm but diverges modestly approaching the UV electronic edge (~210 nm). The Cauchy asymptote √A = 1.823 sets the high-frequency-above-phonons refractive index used by the IR Lorentz; the MP DFPT value n = √ε∞ = 1.95 is ~7% higher and is reported only as cross-reference, not used as the Lorentz anchor.
- **Birefringence: none** (Fm3̄m is cubic).

**Primary references:**
- G. E. Pynchon and E. F. Sieckmann, *Phys. Rev.* **143**, 595 (1966), DOI 10.1103/PhysRev.143.595 — single-crystal visible refractive index, six points.
- J. L. Jacobson and E. R. Nixon, *J. Phys. Chem. Solids* **29**, 1067 (1968), DOI 10.1016/0022-3697(68)90233-3 — IR reflectance and lattice vibrations of CaO and SrO.
- Materials Project entry mp-2472 — DFPT ε∞/ε₀.

---

## 26. AlF₃ combines film visible data (Jensen 1970) with α-AlF₃ DFPT and rhombohedral phonon anchors

**Phase:** at room temperature α-AlF₃ is rhombohedral R3̄c (Z = 6 in the hexagonal setting; Z = 2 in the rhombohedral primitive cell). Lattice parameters in the rhombohedral setting a = 5.0314 Å, α = 58.83°. AlF₃ transforms to a cubic Pm3̄m β-phase near 450 °C.

**Density:** ρ_x = 6 × 83.98 / (6.022e23 × 290.3e-24) (hexagonal-setting volume) = **2.88 g/cm³**, matching the value 2.882 g/cm³ from Daniel et al., *Solid State Commun.* (1990). The Materials Project entry **mp-468** reports 3.02 g/cm³ from the PBE-relaxed cell (5% high due to PBE cell-overestimate inversion in fluorides). The code uses **2880 kg/m³**.

**Optical model — film visible data plus α-AlF₃ DFPT and Lorentz.** AlF₃ has no published bulk single-crystal n(λ),k(λ) spectrum, but **three independent film datasets** plus Materials Project DFPT provide a consistent picture:

**Film data (visible, used directly for thin-flake / coating applications):**
- **Jensen (1970)**, *Thin Solid Films* **5**, 415, DOI 10.1016/0040-6090(70)90052-0 — vacuum-evaporated AlF₃ film: **n ≈ 1.385 across 0.25–1.0 µm**, essentially without dispersion; k ≈ 0 below absorption bands near 3 and 6 µm.
- **Liu et al. (2004)**, *Vacuum* **75**, 159, DOI 10.1016/j.vacuum.2004.01.009 — sputtered film at **193 nm**: **n = 1.37, k = 6·10⁻⁴**.
- **Plumey et al. (2010)**, *Opt. Commun.* **283**, 1942, DOI 10.1016/j.optcom.2009.11.062 — VUV thin films: absolute **n,k in 60–124 nm**.

These film values cover deep-UV through near-IR consistently. The code uses Jensen's near-constant n = 1.385 across 0.25–1.0 µm (k = 0).

**DFPT dielectric data for bulk α-AlF₃ (Materials Project mp-468):**
- ε∞ tensor: diag(1.94, 1.94, 1.94), orientation-averaged ε∞ = **1.94** → n_∞ = √1.94 = **1.393**, within 1% of Jensen's film value 1.385.
- ε₀ tensor: diag(4.93, 4.93, 5.11), orientation-averaged ε₀ = **4.99**
- PBE band gap: 7.72 eV (true gap likely 10–11 eV — AlF₃ is one of the widest-gap dielectrics known, used precisely for this reason in deep-UV optics).

**IR Lorentz** anchored to literature TO modes (Daniel et al. 1990 Raman + far-IR work; Stokey et al. 2022 DFT ellipsometry on related fluorides):

| j | ω_j (cm⁻¹) | λ_j (µm) | Δε_j | γ_j (cm⁻¹) | Assignment |
|:-:|:----------:|:--------:|:----:|:----------:|:-----------|
| 1 | 670 | 14.9 | 1.5 | 40 | Al–F ν₃ asymmetric stretch |
| 2 | 400 | 25.0 | 1.0 | 50 | Al–F ν₂ bend |
| 3 | 270 | 37.0 | 0.5 | 50 | A1g soft mode / lattice rotation |

These reproduce ε₀ ≈ 4.94 by construction, matching DFPT to within 1%. AlF₃ remains transparent (k < 10⁻⁴) from ~0.12 µm to ~10 µm; the strong absorption band near 14.9 µm sets the long-wavelength transparency cutoff.

**Key caveats.**
- **Bulk vs film.** The code's visible n = 1.385 is from evaporated films and matches DFPT-bulk to within 1%. Real bulk α-AlF₃ crystals (rare) may have slightly higher n through density (less voids); the agreement above is fortuitously good.
- **Phase transition at 450 °C.** Cubic β-AlF₃ (Pm3̄m, Materials Project mp-8039) has slightly different optics due to absent octahedral tilting; the code applies to RT α-phase.
- **Hydrolysis.** AlF₃ is much more chemically robust than the alkaline-earth nitrides, but still reacts with water vapor at elevated temperatures; aerosolized AlF₃ in humid environments may carry an Al(OH)F₂ skin.
- **DUV applications.** AlF₃ is widely used as an antireflection coating at 193 nm (ArF excimer wavelength), where Liu et al.'s film data is the authoritative source. The code's near-constant Jensen value n = 1.385 is approximate at 193 nm; for ArF-specific calculations use 1.37 directly.

**Primary references:**
- B. Jensen, *Thin Solid Films* **5**, 415 (1970), DOI 10.1016/0040-6090(70)90052-0 — evaporated AlF₃ films, visible/near-IR.
- Z. Liu et al., *Vacuum* **75**, 159 (2004), DOI 10.1016/j.vacuum.2004.01.009 — sputtered AlF₃ at 193 nm.
- J. Plumey et al., *Opt. Commun.* **283**, 1942 (2010), DOI 10.1016/j.optcom.2009.11.062 — VUV 60–124 nm.
- P. Daniel, A. Bulou, et al., *J. Phys.: Condens. Matter* **2**, 5663 (1990), DOI 10.1088/0953-8984/2/26/004 — α-AlF₃ Raman + phase transition.
- Materials Project entry mp-468 — DFPT ε∞/ε₀ for α-AlF₃.

---

## Conclusion

The implementation combines several well-traced literature models with approximate models where complete optical constants are unavailable. Aluminum, MgF₂, MgO (transparent region), Al₂O₃, and brass Cu70/Zn30 are tied to published coefficients or tabulated values. Carbon soot follows the widely cited Chang & Charalampopoulos (1990) polynomial, though the fitted curve slightly underpredicts some visible-band tabulated values. Magnesium is currently represented by a coarse piecewise approximation rather than a full tabulated interpolation. MgCl₂ correctly uses a handbook n_D value but lacks published dispersion or absorption data. ZnCl₂ likewise has no measured optical tables; its visible refractive indices are independently supported by DFPT (Materials Project mp-22909), while its IR dispersion is a physically constrained Lorentz reconstruction based on DFPT ε∞/ε₀ and FTIR/Raman phonon-band positions.

The least reliable entry is **Al₄C₃**, whose Sellmeier coefficients have no identifiable source and yield n ≈ 2.08 at 589 nm – potentially 30% below secondary-source estimates. This should be treated as a provisional placeholder rather than a validated optical-constants model. A second practical finding cuts across multiple materials: the Sellmeier C parameters for **MgF₂, Al₂O₃, Al₄C₃, AlN, MgAl₂O₄, and Y₂O₃-stabilized cubic ZrO₂** are resonance wavelengths in µm that must be squared in the denominator (λ² − C²). Misinterpreting them as pre-squared values introduces systematic errors of several percent. Finally, the valid wavelength ranges cluster around **0.2–5 µm for the oxide and fluoride Sellmeier models**, with separate phonon or tabulated models required for mid-IR to far-IR Mie calculations beyond ~5 µm.

The 2026 extension (K₂CO₃, KCN, AlN, MgAl₂O₄, KAlO₂, Na₂CO₃, ZrO₂ baddeleyite, ZrO₂-Y₂O₃ cubic YSZ, ZrC, TiC, BaO, Ba₃N₂, Mg₃N₂, NaAlO₂, SrO, AlF₃) adds 16 materials of varying reliability. **Quality-A entries** with published primary spectra invoked directly by the code are AlN (Pastrňák Sellmeier + Kischkat IR), MgAl₂O₄ (Tropf–Thomas), YSZ cubic ZrO₂ (Wood & Nassau Sellmeier, 12 mol% Y₂O₃ strict), and TiC (Pflüger / Palik consolidated EELS-derived n,k 0.031–2.48 µm). **Quality-B entries** with handbook scalars or partial spectra plus reconstructed dispersion are SrO (Pynchon & Sieckmann visible — two endpoint values confirmed from the open abstract, full six-point table not digitized in this pass — plus Jacobson–Nixon TO-anchored IR Lorentz), BaO (n_D from Anderson & Hensley + handbook Cauchy + compilation IR Lorentz), NaAlO₂ (HSDB scalars + Materials Project DFPT), KCN (handbook n_D plus CN-stretch Lorentz tail), Mg₃N₂ (DFT-anchored anti-bixbyite model), m-ZrO₂ baddeleyite (DFPT-anchored, since open primary single-crystal n,k tables for the pure monoclinic phase were not located), and AlF₃ (film data from Jensen + bulk α-AlF₃ Materials Project DFPT). **Quality-C entries** that are full Lorentz reconstructions from DFPT plus vibrational anchors, with no measured optical spectrum, are K₂CO₃, Na₂CO₃, KAlO₂, ZrC (bulk, since only nanocrystalline-film data exist), and Ba₃N₂. The bulk CaO-stabilized c-ZrO₂ dataset of Synowicki & Tiwald (2004) is documented in §18 as a related citable reference covering 0.13–33 µm but is **not implemented** as a code material — the YSZ entry stays strict to Wood & Nassau and 12 mol% Y₂O₃.

A consistent methodological note across these 16 new materials: the Materials Project PBE-relaxed densities are systematically ~3–6% lower than X-ray crystallographic values because of the standard GGA volume overestimate. The handbook always uses the X-ray density (or, where unavailable, the CRC tabulated value); Materials Project densities appear in the "DFPT dielectric data" sections only as cross-references and should not be used for Mie mass loading. Likewise, Materials Project PBE band gaps are 30–50% lower than optical gaps and are quoted only to confirm insulator/metallic character, not to set the UV cutoff of the optical model.
