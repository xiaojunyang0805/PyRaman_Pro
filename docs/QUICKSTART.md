# PyRaman Pro — QuickStart

> PyRaman Pro is a single-window **Spectrum Processor**: a Spectra list on the
> left, the plot in the middle, and a step-by-step Processing panel on the right.

## Running

**From source**
```
cd src
python PyRamanGUI.py
```

**Executable** (after building — see README)
```
dist\PyRaman_Pro\PyRaman_Pro.exe
```
Requirements: Python ≥3.10 with `PyQt5, numpy, scipy, matplotlib, pandas, openpyxl, pybaselines` (`pip install -r requirements.txt`).

---

## The workflow: Open → process → export

### 1. Open spectra
Click **➕ Open files…** (or **File ▸ Open files…**, Ctrl+O) and pick one or more
`.csv` / `.txt` / `.dat` / `.xlsx` files. Each file's first numeric column is the
Raman shift (X); every other numeric column becomes a spectrum (Y). Delimiter,
decimal comma, BOM/encoding and a header row are auto-detected.

Loaded spectra appear in the **Spectra** list (left). Click one to make it
active; Ctrl/Shift-click several to **overlay** them for comparison.

### 2. Auto Analyze (optional first pass)
**Auto Analyze** picks a baseline and detects peaks automatically. Good for
smooth fluorescence; for strongly structured backgrounds (see polystyrene below)
the manual Strength slider works better.

### 3. Tune the pipeline (right panel)
Each step is a numbered on/off card with a live preview:

| Step | Does | Tips |
|---|---|---|
| ① Range (crop) | limit the wavenumber window | e.g. start at 200 cm⁻¹ |
| ② Despike | remove single-point cosmic spikes | raise/lower σ |
| ③ Smooth | Savitzky–Golay | window 9, poly 3 is a good start |
| ④ Baseline | subtract background | **arPLS** (robust) or **airPLS** (flexible); **Strength** = λ |
| ⑤ Normalize | scale for comparison | Min-Max / Max=1 / Area=1 |
| ⑥ Peaks | detect bands | **% of strongest** — higher = fewer, cleaner |

**Strength (λ):** the slider spans λ=10¹…10⁹ (the live value shows next to it).
*Lower λ = more flexible baseline* that follows broad humps; *higher λ = stiffer*.

**IS-Score** (Baseline step) is an objective 0–1 quality score for the current
baseline (green ≥0.8 / amber ≥0.6 / red below). **Optimize (IS-Score)** searches
λ to maximise it. See `doc/Polystyrene_demo.md` for what IS-Score is and its one
caveat.

### 4. Manual anchor baseline (OriginLab-style)
In **Baseline** click **✎ Edit anchors**:
- Anchors auto-seed on the **low-curvature background points** (like OriginLab's
  2nd-derivative method) and a **B-spline** is drawn through them.
- **Anchor count** (default 40), **Flatness %** (curvature threshold), and
  **Find anchors** re-seed; **Interpolation** = B-spline / PCHIP / Linear.
- **Drag** a handle up/down (x is locked by default) · **double-click** to add ·
  **right-click** to delete · **Snap to valley** keeps the baseline under the signal.

### 5. Navigate the plot
- **Mouse wheel** = zoom (centred on cursor); **Shift + wheel** = zoom Y.
- **Left-drag** = pan. **Reset view** (toolbar) = zoom back out.
- Zoom in to place/drag anchors precisely.

### 6. Fit peaks, preview, export
- **Tools ▸ Fit peaks** fits overlapping bands jointly (Lorentzian) → Position,
  Height, FWHM, Area in the table.
- **Preview Result** = clean full-scale view of the corrected spectrum.
- **Export Data** (CSV x,y) / **Export Figure** (PNG/PDF/SVG).

---

## Polystyrene example (file `examples/data/polystyrene.csv`)
Recommended: **Smooth (9, 3) + airPLS, Strength λ≈100** (slide Strength low until
the peaks emerge from the humps). Recovers the polystyrene fingerprint —
**~1001 (strongest), 1031, 620, 795, 1331, 1450, 1602, 2852, 2904, 3054 cm⁻¹** —
flat between 1700–2800. Full walkthrough: `doc/Polystyrene_demo.md`.

## Shortcuts
| Key | Action |
|---|---|
| Ctrl+O | Open files |
| Ctrl+Q | Exit |
| Mouse wheel | Zoom (Shift = Y) |
| Left-drag | Pan |
