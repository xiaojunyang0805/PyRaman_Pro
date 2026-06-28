# PyRamanGUI vs OriginLab - Feature Comparison

> **Update 2026-06-28 — single-window app, anchor parity reached.**
> PyRamanGUI is now a single-window **Spectrum Processor** (Spectra list · plot ·
> Processing panel). The OriginLab-style shell (folders, spreadsheets, plot/text/
> mapping windows, `.rmn` projects) was **removed** — ignore the "Spreadsheet /
> Plot window / folder tree" descriptions below; they are historical.
> - **Anchor editing now mirrors OriginLab directly:** *✎ Edit anchors* seeds
>   anchors on **low-curvature (2nd-derivative) background points** (OriginLab's
>   "1st & 2nd Derivative (zeroes)"), interpolates with a **B-spline** (default;
>   PCHIP/Linear also available), with **Anchor count**, **Flatness %**, **Find
>   anchors**, drag/add/delete, and **mouse-wheel zoom** to edit in detail.
> - **Baseline methods** trimmed to the two that matter for Raman: **arPLS**
>   (robust) and **airPLS** (flexible), with a Strength (λ) slider spanning 10¹–10⁹.
> - **Objective IS-Score** (ACS Omega 2026) rates baseline quality live — a
>   capability OriginLab does not have. See `doc/Polystyrene_demo.md`.
> - **Peak fitting** reports FWHM + area with joint fitting of overlapping bands.
> - Workflow is now **➕ Open files… → tune the pipeline → Export**. See
>   `PyRamanGUI/QUICKSTART.md`.

## Your OriginLab Workflow (from doc/241217)

Based on your documented workflow for processing Raman spectra in OriginLab:

### Workflow Steps in OriginLab:
1. **Plot raw Raman data**
2. **Smooth the data**
3. **Baseline correction** using Peak Analyzer:
   - Goal: Subtract Baseline
   - Method: 1st and 2nd Derivative (zeroes)
   - Selecting baseline anchor points (User Defined)
   - Interpolation: BSpline
   - Manual adjustment of anchor points (Add/Modify/Delete)
   - Subtract baseline
4. **Plot subtracted data** in new workspace
5. **Compare** with reference (e.g., polystyrene reference)

---

## Can PyRamanGUI Do This? ✅ YES!

### Feature-by-Feature Comparison

| OriginLab Workflow Step | PyRamanGUI Equivalent | Status | Notes |
|------------------------|----------------------|--------|-------|
| **1. Plot raw data** | Spreadsheet → Plot Window | ✅ **YES** | Import data, select columns, create plot |
| **2. Smoothing** | Analysis → Smoothing | ✅ **YES** | Multiple methods available |
| **3. Baseline correction** | Analysis → Baseline Correction | ✅ **YES** | 10+ methods including splines |
| **Anchor points (manual)** | Spectrum Processor → ✎ Edit anchors | ✅ **YES** | Draggable points, add/delete, snap-to-valley, live preview (PCHIP) |
| **BSpline interpolation** | Available in spline methods | ✅ **YES** | Univariate Spline, GCV Spline |
| **Plot subtracted data** | Automatic in Plot Window | ✅ **YES** | Shows both original and corrected |
| **Compare to reference** | Database for measurements | ✅ **YES** | Store reference peak positions |

---

## Detailed Comparison

### 1. Data Import & Plotting ✅

**OriginLab**: Import → Plot
**PyRamanGUI**:
- Create Spreadsheet (File → New → Spreadsheet)
- Load Data (File → Load Data)
- Select X and Y columns
- Right-click column → "Plot in new Plot Window"

**Verdict**: ✅ **Equivalent functionality**

---

### 2. Smoothing ✅

**OriginLab**: Analysis → Signal Processing → Smoothing
**PyRamanGUI**: Analysis → Smoothing

**Available Methods in PyRamanGUI**:
- **Whittaker smoothing**
- **Spline smoothing**
- **Window methods** (Moving average, Savitzky-Golay)

**Verdict**: ✅ **Similar or better** - Multiple sophisticated methods

---

### 3. Baseline Correction ✅

This is where PyRamanGUI really shines!

**OriginLab**:
- Peak Analyzer → Subtract Baseline
- Methods: User-defined anchor points, derivatives
- Interpolation: BSpline, Linear, etc.

**PyRamanGUI**: Analysis → Baseline Correction

**Available Methods** (10+ algorithms):

#### **A. Whittaker-based Methods** (Recommended)
1. **Asymmetric Least Square (ALS)** ⭐
   - Most popular method
   - Automatic, no anchor points needed
   - Parameters: `p` (asymmetry), `lambda` (smoothness)

2. **Improved ALS (imALS)**
   - Enhanced version of ALS

3. **airPLS** (Adaptive Iteratively Reweighted)
   - Fully automatic
   - Good for complex baselines

4. **arPLS** (Asymmetrically Reweighted)
   - Robust to peak interference

5. **drPLS** (Doubly Reweighted)
   - High quality results

6. **derpsALS** (Derivative Peak-Screening)
   - Uses derivatives (similar to OriginLab!)

#### **B. Spline Methods** ⭐ (Like OriginLab BSpline)
1. **Univariate Spline**
   - User-defined anchor points (ROI)
   - Direct equivalent to OriginLab BSpline!
   - Parameters: `s` (smoothing), `roi` (anchor regions)

2. **GCV Spline**
   - Generalized Cross-Validation
   - User-defined anchor points

#### **C. Polynomial Methods**
1. **Polynomial** (with ROI)
   - User-defined baseline regions
   - Fit polynomial through anchor points

2. **Polynomial** (without ROI)
   - Automatic polynomial fitting

#### **D. Miscellaneous**
1. **Rubberband**
   - Convex hull method
   - Fully automatic

2. **Rolling Ball**
   - Morphological method
   - Parameter: half window size

**Verdict**: ✅ **SUPERIOR** - More methods, more flexibility than OriginLab!

---

### 4. Interactive Anchor Point Selection ⚠️

**OriginLab**:
- Visual selection of anchor points
- Add/Modify/Delete interactively
- See preview of baseline

**PyRamanGUI**:
- **Spline methods**: Define ROI (regions of interest) as anchor areas
- **Most methods**: Automatic (no anchor points needed)
- Less interactive than OriginLab for manual selection

**Verdict**: ⚠️ **Different approach**
- OriginLab: More interactive point-by-point control
- PyRamanGUI: More algorithmic, less manual intervention
- **Recommendation**: Use Univariate Spline with ROI for closest match

---

### 5. Reference Comparison ✅

**OriginLab**: Manual overlay of reference spectra

**PyRamanGUI**:
- Tools → Database for measurements
- Add reference peak positions
- Mark reference lines on plots
- Overlay multiple spectra in same plot window

**Verdict**: ✅ **Equivalent** - Can overlay and compare

---

## PyRamanGUI Advantages Over OriginLab

### 1. **More Baseline Correction Methods**
- OriginLab: ~5 methods
- PyRamanGUI: 10+ advanced methods
- State-of-the-art algorithms (ALS, airPLS, etc.)

### 2. **Better for Batch Processing**
- Create analysis routines (drag-and-drop)
- Apply to multiple spectra at once
- Reproducible workflows

### 3. **Project Organization**
- Folder structure for organizing data
- Save entire project as single .rmn file
- Text windows for notes

### 4. **Free and Open Source**
- No license costs
- Modify code if needed
- Community contributions

### 5. **Multivariate Analysis**
- PCA (Principal Component Analysis)
- NMF (Non-negative Matrix Factorization)
- Good for complex datasets

### 6. **Peak Fitting**
- Multiple line shapes:
  - Lorentzian
  - Gaussian
  - Voigt
  - Pseudo-Voigt
  - Breit-Wigner-Fano
- Batch fitting across spectra

---

## OriginLab Advantages Over PyRamanGUI

### 1. **More Interactive Baseline Selection**
- Point-by-point anchor selection
- Real-time visual feedback
- Easier for beginners

### 2. **More Mature UI/UX**
- Polished interface
- Better documentation
- More intuitive for some users

### 3. **Broader Functionality**
- Not just Raman - all kinds of data
- More general-purpose tool
- More analysis types

### 4. **Professional Support**
- Official customer support
- Regular updates
- Training materials

---

## Recommended PyRamanGUI Workflow (Matching Your OriginLab Process)

### For Your Polystyrene Analysis:

**Step 1: Import Raw Data**
```
File → New → Spreadsheet
File → Load Data
Select your Raman .txt/.csv file
Set delimiter (Tab/Space/Comma)
```

**Step 2: Create Initial Plot**
```
Select Y column (Intensity)
Right-click → "Plot in new Plot Window"
```

**Step 3: Smooth Data**
```
In Plot Window:
Analysis → Smoothing
Choose method (e.g., Savitzky-Golay)
Adjust parameters
Apply
```

**Step 4: Baseline Correction**

**Option A: Closest to Your OriginLab Method**
```
Analysis → Baseline Correction
Select "Univariate Spline"
Parameters:
  - s: smoothing parameter (start with 1e0)
  - roi: define baseline regions, e.g., [[200,300], [3500,4000]]
Apply
```

**Option B: Automatic (Recommended)**
```
Analysis → Baseline Correction
Select "Asymmetric Least Square (ALS)"
Parameters:
  - p: 0.001 (default)
  - lambda: 10000000 (default)
Apply
(No anchor points needed - fully automatic!)
```

**Step 5: View Results**
```
Baseline-corrected spectrum appears in plot
Export plot: Toolbar → Save icon
Export data: Spreadsheet → File → Export
```

**Step 6: Compare with Reference**
```
Tools → Database for measurements
Add polystyrene reference peaks
Or: Import reference spectrum and overlay
```

---

## Key Differences in Philosophy

### OriginLab Approach:
- **Manual/Interactive**: User controls anchor points
- **Visual/Intuitive**: See and adjust in real-time
- **General-purpose**: Works for any data type

### PyRamanGUI Approach:
- **Algorithmic/Automatic**: Advanced algorithms do the work
- **Less manual intervention**: More reproducible
- **Raman-specific**: Optimized for spectroscopy

---

## Practical Recommendations

### For Your Workflow:

1. **Try ALS (Asymmetric Least Square) first** ⭐
   - Most popular baseline correction method
   - Fully automatic
   - Usually gives excellent results
   - Adjust `lambda` (higher = smoother baseline)

2. **If you need manual control like OriginLab:**
   - Use **Univariate Spline** with ROI
   - Define regions where baseline should be (no peaks)
   - Similar to your anchor point method

3. **For polystyrene analysis:**
   - Polystyrene has well-defined peaks
   - ALS or airPLS should work perfectly
   - Less need for manual anchor points

### Suggested Testing Workflow:

```
1. Import your polystyrene data
2. Create plot
3. Try baseline correction with ALS (default parameters)
4. Compare result with your OriginLab output
5. If needed, adjust parameters or try other methods
6. Once satisfied, save as analysis routine for reuse
```

---

## Feature Support Summary

| Feature Category | OriginLab | PyRamanGUI | Winner |
|-----------------|-----------|------------|--------|
| Data Import | ✅ | ✅ | Tie |
| Plotting | ✅ | ✅ | Tie |
| Smoothing | ✅ | ✅ | Tie |
| Baseline Correction | ✅ | ✅✅ | **PyRamanGUI** (more methods + per-spectrum λ search) |
| Interactive Anchors | ✅✅ | ✅✅ | **Tie** (draggable anchors + snap-to-valley, since 2026-06-20) |
| Peak Fitting | ✅ | ✅ | Tie (FWHM/area, joint fit of overlaps) |
| Batch Processing | ✅ | ✅✅ | **PyRamanGUI** (better routines) |
| Multivariate Analysis | ⚠️ | ✅ | **PyRamanGUI** (PCA/NMF) |
| UI/UX Polish | ✅✅ | ✅ | **OriginLab** (more mature) |
| Cost | ❌ $$$$ | ✅ FREE | **PyRamanGUI** |

---

## Bottom Line

### Can PyRamanGUI replace OriginLab for your Raman workflow?

# ✅ **YES!**

PyRamanGUI can do everything in your documented OriginLab workflow, often with **more sophisticated algorithms**.

### Key Points:

1. **Baseline correction**: PyRamanGUI has MORE methods than OriginLab
2. **Smoothing**: Equivalent functionality
3. **Plotting**: Similar capabilities
4. **Interactive anchors**: Different approach (more automatic, less manual)
5. **Overall**: PyRamanGUI is excellent for Raman analysis, especially with modern algorithms

### Recommendation:

✅ **Test PyRamanGUI with your polystyrene data**
✅ **Start with ALS baseline correction** (easier than manual anchors)
✅ **You'll likely get equal or better results with less manual work**

---

## Next Steps for You

1. **Tomorrow**: Import your polystyrene Raman data into PyRamanGUI
2. **Try the workflow** I outlined above
3. **Compare results** with your OriginLab output
4. **Document**: Which methods work best for your data
5. **Create routine**: Save successful workflow for reuse

Let me know how it goes! 🎉

---

**Last Updated**: 2026-06-20
