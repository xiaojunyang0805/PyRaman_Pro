# Polystyrene demo — baseline correction, peaks, and the IS-Score

A worked example using **`examples/data/polystyrene.csv`** ("Polystyrene beads on
glass slide"), showing how to recover the polystyrene Raman fingerprint in the
PyRaman Pro **Spectrum Processor** and what each method does. The target is the
published polystyrene Raman reference (ASTM E1840; see any standard polystyrene
Raman spectrum).

---

## 1. What the raw data looks like

The genuine polystyrene peaks are **sharp but small**, sitting on a **large,
wavy/structured background** (period ~600–700 cm⁻¹) from glass fluorescence and
etaloning — *not* a smooth fluorescence slope. So the baseline must be **flexible
enough to follow the humps** while leaving the narrow peaks intact.

Reference polystyrene bands (ASTM E1840):

| cm⁻¹ | assignment | cm⁻¹ | assignment |
|---|---|---|---|
| 620 | ring deformation | 1450 | CH₂ scissor |
| 795 | CH₂ rock | 1583 | ring |
| **1001** | **ring breathing (strongest)** | **1602** | ring stretch |
| 1031 | ring C–H in-plane | 2852 | CH₂ sym |
| 1155 / 1183 | C–C / ring C–H | 2904 | C–H |
| 1331 | CH deformation | **3054** | aromatic C–H |

Between **1700–2800 cm⁻¹** the corrected spectrum should be **flat at zero**.

---

## 2. Recommended procedure (automatic)

1. **➕ Open files…** → select `examples/data/polystyrene.csv`.
2. Enable **③ Smooth** — window **9**, poly **3**.
3. Enable **④ Baseline** — Method **airPLS (flexible)** (or arPLS).
4. Drag **Strength** *down* to about **λ = 1×10²** (slider near the low end).
   Watch the broad humps flatten and the sharp peaks (1001, 1602, 3054) rise out
   of the background. The live λ value shows next to the slider.
5. Enable **⑥ Peaks** to mark the bands; **Tools ▸ Fit peaks** for FWHM/area.
6. **Preview Result** → **Export Data / Figure**.

This recovers ~10 of the 13 reference bands with a flat baseline, dominated by
1001 / 1602 / 3054 — matching the reference.

> **Why low λ?** The structured glass background varies on a ~600 cm⁻¹ scale.
> A stiff baseline (high λ, e.g. 10⁴) arches over the humps and leaves them in
> the result; λ≈10² lets the baseline follow the humps. The Strength slider now
> reaches this low range (it used to stop at 10⁴).

---

## 3. Recommended procedure (manual anchors, OriginLab-style)

Use this when you want hands-on control, or to mirror the OriginLab Peak
Analyzer workflow.

1. Smooth as above; enable **④ Baseline**.
2. Click **✎ Edit anchors**. Anchors auto-seed on the **low-curvature
   (background) points** — the same idea as OriginLab's *"1st & 2nd Derivative
   (zeroes)"* method — and a **B-spline** is drawn through them.
3. Tune density with **Anchor count** (default 40) and **Flatness %** (curvature
   threshold; lower = stricter = fewer, flatter anchors). **Find anchors**
   re-seeds.
4. **Zoom in** (mouse wheel) on a region and **drag** individual anchors up/down
   to fine-tune; **double-click** to add, **right-click** to delete. **Reset
   view** zooms back out.
5. Choose **Interpolation: B-spline** (smooth, like OriginLab) / PCHIP / Linear.

The dense anchors + B-spline follow the broad background and skip the sharp
peaks, giving the same fingerprint as the automatic route.

---

## 4. The IS-Score (baseline-quality metric)

The **IS-Score (Integrity Spectrum Score)** shown in the Baseline step is a
**clean-room re-implementation** of:

> S. Innocente, A. Visentin, S. Andersson-Engels, K. Komolibus, R. Gautam,
> *Automated Baseline Correction Evaluation Score for Raman Spectroscopy*,
> ACS Omega 2026, 11, 25057–25068. DOI: 10.1021/acsomega.5c09870

**What it is.** A single number in **[0, 1]** (higher = better) that judges a
baseline correction **without any ground truth** — from just the raw spectrum,
the corrected spectrum, and the axis. It encodes what an expert checks by eye.

**How it works.** `IS-Score = 1 − Σ(penalties)` over five blocks:

1. **Single peak/dip** — the baseline must not cut into peaks (overfit) or sit
   off the dips (underfit), judged against each band's prominence.
2. **Band region (RRP)** — extends that across each peak's full region via a
   baseline-scaled "Raman Region Prominence".
3. **AUC** — penalises global over-correction (negative area where the baseline
   rose above the signal).
4. **Mean ratio** — penalises dips left **above** the baseline (leftover
   background = underfitting).
5. **Intensity** — penalises points where the baseline exceeds the spectrum
   (beyond noise).

It is computed live on every change; colour: green ≥0.8, amber ≥0.6, red below.
**Optimize (IS-Score)** searches λ to maximise it.

**Important caveat.** IS-Score was validated on *smooth* fluorescence (blood
cells). On this polystyrene data's **structured/wavy** background it slightly
prefers the *wrong* stiff baseline (λ≈10⁴, IS≈0.76, which leaves the humps) over
the correct flexible one (λ≈10², IS≈0.75, which recovers the peaks). So for
structured backgrounds, **trust the visible peaks and the manual Strength slider
rather than maximising IS-Score**. For ordinary smooth fluorescence, IS-Score and
Optimize work well.

---

## 5. How this maps to the OriginLab method

| OriginLab Peak Analyzer | PyRamanGUI Spectrum Processor |
|---|---|
| Smooth | ③ Smooth (Savitzky–Golay) |
| Anchor finding: "1st & 2nd Derivative (zeroes)", ~60 points | ✎ Edit anchors → low-curvature seeding, Anchor count + Flatness |
| Interpolation: BSpline | Interpolation: **B-spline** (default) |
| Add / Modify / Delete anchors | drag / double-click / right-click (+ zoom) |
| Subtract baseline | live, automatic |
| (none) | objective **IS-Score** + automatic arPLS/airPLS |

Both routes produce the same clean polystyrene spectrum; PyRamanGUI adds the
automatic baseline + objective scoring on top.
