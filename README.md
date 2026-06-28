# PyRaman Pro

A simple, spectrum-centric desktop app for **Raman spectrum processing** —
baseline correction, smoothing, peak detection/fitting — with an objective,
ground-truth-free baseline-quality metric (**IS-Score**).

It opens straight into a single window: a **Spectra** list on the left, the plot
in the middle, and a step-by-step **Processing** panel on the right.

```
Open files → Range → Despike → Smooth → Baseline → Normalize → Peaks → Export
```

## Download (Windows)

Grab the latest ready-to-run build from the
**[Releases page](https://github.com/xiaojunyang0805/PyRaman_Pro/releases/latest)**:

1. Download `PyRaman_Pro-windows-x64.zip`.
2. Extract it anywhere.
3. Run **`PyRaman_Pro\PyRaman_Pro.exe`** — no Python install needed.

(No release yet? Build it yourself — see [Build a Windows executable](#build-a-windows-executable).)

## Features

- **Open** one or more `.csv` / `.txt` / `.dat` / `.xlsx` spectra (auto delimiter,
  decimal-comma, encoding, header). Overlay several to compare.
- **Baseline correction** — **arPLS** (robust) and **airPLS** (flexible), with a
  Strength (λ) slider spanning 10¹–10⁹, plus a manual **anchor + B-spline** mode
  (OriginLab-style: low-curvature/derivative-zero seeding, Find anchors, Flatness,
  drag/add/delete, mouse-wheel zoom).
- **IS-Score** — a live 0–1 baseline-quality score, a clean-room re-implementation
  of Innocente et al., *Automated Baseline Correction Evaluation Score for Raman
  Spectroscopy*, ACS Omega 2026, 11, 25057 (DOI 10.1021/acsomega.5c09870), with a
  per-spectrum λ optimizer.
- **Peaks** — noise-adaptive detection + Lorentzian fitting (FWHM, area; joint fit
  of overlapping bands).
- **Plot** — mouse-wheel zoom (Shift = Y), drag pan, Reset view.
- **Export** data (CSV) and figures (PNG/PDF/SVG).

## Install & run

Requires Python ≥ 3.10.

```bash
pip install -r requirements.txt
cd src
python PyRamanGUI.py
```

### Build a Windows executable
```bash
python -m PyInstaller PyRaman_Pro.spec --noconfirm --clean
# -> dist/PyRaman_Pro/PyRaman_Pro.exe
```

## Try it

```bash
cd src
python PyRamanGUI.py
# File ▸ Open files… → examples/data/polystyrene.csv
```
Then follow **[docs/Polystyrene_demo.md](docs/Polystyrene_demo.md)** to recover
the polystyrene fingerprint (Smooth 9/3 + airPLS λ≈100).

## Documentation

- [QuickStart](docs/QUICKSTART.md) — the full workflow.
- [Polystyrene demo](docs/Polystyrene_demo.md) — worked example + IS-Score explained.
- [OriginLab comparison](docs/OriginLab_comparison.md).
- [autoRaman package](src/autoRaman/README.md) — the GUI-free analysis core + IS-Score.

## Tests

```bash
cd src
python -m pytest autoRaman/tests/ -q
```

## Credits & License

Derived from **PyRamanGUI** by Simon Brehm
([gitlab.com/brehmsi/PyRamanGUI](https://gitlab.com/brehmsi/PyRamanGUI)), reworked
into the single-window Spectrum Processor with the IS-Score and OriginLab-style
anchor editing. Licensed under the **Apache License 2.0** — see [LICENSE](LICENSE).

The IS-Score is an independent re-implementation of the *method* described in the
ACS Omega 2026 paper (no upstream source code was used); please cite that paper if
you use the score.
