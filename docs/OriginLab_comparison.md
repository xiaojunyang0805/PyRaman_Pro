# PyRaman Pro vs OriginLab — baseline correction workflow

This page maps a typical OriginLab workflow for Raman baseline correction onto
PyRaman Pro, and is honest about where each tool is stronger. The worked
example on the same kind of data is in [Polystyrene_demo.md](Polystyrene_demo.md).

## A typical OriginLab workflow

A common way to process a Raman spectrum in OriginLab (Peak Analyzer):

1. **Plot** the raw spectrum.
2. **Smooth** it (Analysis ▸ Signal Processing ▸ Smooth).
3. **Subtract the baseline** with Analysis ▸ Peaks and Baseline ▸ Peak Analyzer:
   - Goal: *Subtract Baseline*
   - Baseline mode *User Defined*, anchor finding *1st and 2nd Derivative (zeroes)*
   - Interpolation *BSpline*
   - Adjust anchors by hand (Add / Modify / Delete), then *Subtract*
4. **Plot** the subtracted spectrum.
5. **Compare** it with a reference spectrum (e.g. polystyrene).

## The same steps in PyRaman Pro

| OriginLab step | PyRaman Pro |
|---|---|
| Plot raw data | ➕ **Open files…** — the spectrum is plotted at once (csv/txt/dat/asc/xlsx; delimiter, decimal comma and header auto-detected) |
| Smooth | ③ **Smooth** — Savitzky–Golay (window, polynomial order) |
| Anchors: 1st & 2nd derivative (zeroes) | ④ Baseline ▸ **✎ Edit anchors** — anchors seeded on low-curvature (background) points; **Anchor count**, **Flatness %**, **Find anchors** |
| Interpolation: BSpline | **Interpolation: B-spline** (default), PCHIP or Linear |
| Add / Modify / Delete anchors | double-click to add, drag to move, right-click to delete; mouse-wheel zoom for detail |
| Subtract baseline | live — the corrected spectrum updates on every change |
| Plot subtracted data | shown directly; **Preview Result** for a clean full-scale view |
| Compare with a reference | overlay several spectra (Ctrl/Shift-click in the list) |

## What PyRaman Pro adds

- **Automatic baselines** — arPLS (robust) and airPLS (flexible) with one
  Strength (λ) slider from 10¹ to 10⁹; **Auto Analyze** picks the method and λ.
  On a strongly structured background (e.g. polystyrene on glass) a low λ (≈10²)
  is needed to follow the background humps.
- **IS-Score** — a live 0–1 baseline-quality score that needs no reference
  spectrum, re-implemented from Innocente et al., *ACS Omega* 2026, 11, 25057
  (doi:10.1021/acsomega.5c09870). Treat it as a guide: on structured
  backgrounds its maximum can favour a baseline that leaves broad humps in
  (see the demo).
- **Peak detection and fitting** — noise-adaptive detection and Lorentzian fits
  (position, height, FWHM, area), with overlapping bands fitted jointly.
- **Free and open source** (Apache-2.0); a ready-to-run Windows build.

## Where OriginLab is stronger

- **Batch processing** — Peak Analyzer themes can be applied to many spectra;
  PyRaman Pro processes one spectrum at a time (batch processing is planned).
- **Peak shapes** — OriginLab fits many profiles (Gaussian, Voigt, …);
  PyRaman Pro fits Lorentzians only.
- **Smoothing and baseline options** — OriginLab offers more smoothing filters
  and anchor-finding modes; PyRaman Pro keeps to Savitzky–Golay and the two
  penalised-least-squares methods that suit Raman spectra.
- **General-purpose analysis and graphing**, a mature interface, documentation
  and commercial support.

## Summary

| | OriginLab | PyRaman Pro |
|---|---|---|
| Manual anchor baseline (B-spline) | ✅ | ✅ |
| Automatic baseline (arPLS / airPLS) | — | ✅ |
| Objective baseline-quality score | — | ✅ IS-Score |
| Peak fitting | ✅ many profiles | ✅ Lorentzian |
| Batch processing | ✅ | planned |
| General data analysis / graphing | ✅ | — (Raman-focused) |
| Cost | commercial licence | free, open source |

## References

- airPLS: Z.-M. Zhang, S. Chen, Y.-Z. Liang, *Analyst* 2010, 135, 1138.
- arPLS: S.-J. Baek, A. Park, Y.-J. Ahn, J. Choo, *Analyst* 2015, 140, 250.
- IS-Score: S. Innocente et al., *ACS Omega* 2026, 11, 25057.
- PyRamanGUI (the project PyRaman Pro is derived from): S. Brehm, C. Himcinschi,
  J. Kraus, J. Kortus, *SoftwareX* 2023, 23, 101486, doi:10.1016/j.softx.2023.101486.
