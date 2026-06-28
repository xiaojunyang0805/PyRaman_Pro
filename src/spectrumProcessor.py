"""
spectrumProcessor.py -- spectrum-centric guided processing view.

The friendly front door for Raman analysis: import a file and you immediately
see the plotted spectrum next to a step-by-step Processing panel that mirrors
the canonical Raman pipeline:

    Range  ->  Despike  ->  Smooth  ->  Baseline  ->  Normalize  ->  Peaks

Every step has a sensible default and a plain-language hint. Each control change
updates a live preview on the plot. "Auto Analyze" runs the whole pipeline with
AutoRaman's automatic baseline selection + peak detection; the user can then
tweak any step (semi-automatic).

Automation is additive: this view never removes the spreadsheet / manual
workflow -- it sits on top of the same imported data and reuses the existing
analysis algorithms (autoRaman + scipy).
"""

import os

import numpy as np
from scipy.signal import savgol_filter, find_peaks, medfilt, peak_widths
from scipy.optimize import curve_fit

from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

# Reuse the AutoRaman algorithm core (GUI-free). Guarded so the view still
# works (manual steps only) if the optional package is unavailable.
try:
    from autoRaman import AutoAnalyzer, AutoPeaks, is_score
    from autoRaman.auto_baseline import (
        AutoBaseline, ParameterOptimizer, _arpls, _airpls)
    _AUTORAMAN = True
except Exception as _e:  # noqa: BLE001
    AutoAnalyzer = AutoPeaks = AutoBaseline = ParameterOptimizer = is_score = None
    _AUTORAMAN = False
    print("AutoRaman unavailable in SpectrumProcessor:", _e)


# --- baseline method registry (method label -> callable + which knob it takes)
# The callables return (y_corrected, baseline) and come straight from autoRaman
# so behaviour matches the rest of the app.
if _AUTORAMAN:
    # Two penalized-least-squares methods cover Raman baselines well: arPLS is
    # the robust default; airPLS is more flexible and often best on broad,
    # structured backgrounds (e.g. polystyrene on glass). Both take lambda, so
    # the Strength slider always applies. Polynomial/rolling-ball were dropped.
    BASELINE_METHODS = {
        "arPLS  (recommended)": dict(fn=_arpls, knob="lam"),
        "airPLS  (flexible)": dict(fn=_airpls, knob="lam"),
    }
else:
    BASELINE_METHODS = {}

# knob -> (slider_min, slider_max, default, label, value_from_slider)
# lam is on a log scale (10**n); the others are linear.
KNOB_SPEC = {
    # slider 10..90 -> lambda 1e1..1e9. The LOW end matters: spectra with a
    # broad, structured/wavy background (e.g. glass fluorescence under polystyrene
    # beads) need lambda ~1e2 for the baseline to follow the humps and reveal the
    # sharp Raman peaks. The old floor of 1e4 could not reach that regime.
    "lam": (10, 90, 30, "Strength (λ = 10^n)", lambda s: 10.0 ** (s / 10.0)),
    "poly_order": (1, 9, 3, "Polynomial order", lambda s: int(s)),
    "half_window": (10, 250, 50, "Ball radius", lambda s: int(s)),
}


def despike(y, sigma=5.0, window=5):
    """Remove single-bin cosmic spikes via a median-difference test."""
    if len(y) < window:
        return y
    med = medfilt(y, window if window % 2 else window + 1)
    diff = y - med
    s = np.std(diff)
    if s == 0:
        return y
    out = y.copy()
    spikes = np.abs(diff) > sigma * s
    out[spikes] = med[spikes]
    return out


def _lorentzian(x, xc, h, fwhm):
    """Lorentzian peak (same definition as peakFitting.LorentzFct); area = pi*h*fwhm/2."""
    return h / (1.0 + (2.0 * (x - xc) / fwhm) ** 2)


def _multi_lorentzian(x, *p):
    """Sum of N Lorentzians (p in triples xc,h,fwhm) + linear background (last 2)."""
    a, b = p[-2], p[-1]
    y = a * x + b
    for i in range(0, len(p) - 2, 3):
        y = y + _lorentzian(x, p[i], p[i + 1], p[i + 2])
    return y


def normalize(y, mode):
    if mode.startswith("Min"):
        rng = y.max() - y.min()
        return (y - y.min()) / rng if rng else y
    if mode.startswith("Max"):
        m = y.max()
        return y / m if m else y
    # Area
    a = _trapz(np.abs(y))
    return y / a if a else y


# np.trapz was renamed to np.trapezoid in NumPy 2.0 (trapz removed in 2.x).
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# --- standalone file loader -------------------------------------------------
# Replaces the old spreadsheet import path: parse a data file straight into
# (x, [(name, y), ...]). The first numeric column is the Raman shift (X); every
# remaining numeric column is one intensity spectrum (Y). Robust to delimiter,
# decimal comma, BOM/encoding and header rows -- the same handling the deleted
# DataImportDialog used.

def _read_lines(path):
    """Read text lines robustly across common encodings (handles BOM)."""
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            with open(path, "r", encoding=enc) as fh:
                return fh.read().splitlines()
        except (UnicodeDecodeError, UnicodeError):
            continue
    with open(path, "r", errors="replace") as fh:
        return fh.read().splitlines()


def _sniff_format(lines):
    """Guess (delimiter, decimal) from a sample of data lines."""
    sample = [ln for ln in lines if ln.strip()][:25]
    if not sample:
        return None, "."
    delimiter = None
    for d in (";", "\t"):
        counts = [s.count(d) for s in sample]
        if min(counts) >= 1 and len(set(counts)) == 1:
            delimiter = d
            break
    if delimiter is None:
        counts = [s.count(",") for s in sample]
        if min(counts) >= 1 and len(set(counts)) == 1:
            delimiter = ","
    if delimiter in (";", "\t", None):
        decimal = "," if any("," in s for s in sample) else "."
    else:  # delimiter == ","
        decimal = "."
    return delimiter, decimal


def _looks_like_header(line, delimiter):
    """True if the first line has non-numeric tokens (a column-name row)."""
    toks = line.split(delimiter) if delimiter else line.split()
    n_bad = 0
    for t in toks:
        t = t.strip().replace(",", ".")
        if t == "":
            continue
        try:
            float(t)
        except ValueError:
            n_bad += 1
    return n_bad > 0


def load_spectra_file(path):
    """Load a data file -> (x, [(name, y), ...]).

    Supports .csv/.txt/.dat (auto delimiter + decimal) and .xlsx/.xls.
    Raises ValueError with a readable message if nothing numeric is found.
    """
    ext = os.path.splitext(path)[1].lower()
    base = os.path.splitext(os.path.basename(path))[0]

    if ext in (".xlsx", ".xls"):
        import pandas as pd
        df = pd.read_excel(path)
        num = df.select_dtypes("number")
        if num.shape[1] < 2:
            raise ValueError("Need at least two numeric columns (X and Y).")
        arr = num.to_numpy(dtype=float)
        names = list(num.columns[1:])
    else:
        lines = _read_lines(path)
        delim, decimal = _sniff_format(lines)
        header_names = None
        data_lines = [ln for ln in lines if ln.strip()]
        if data_lines and _looks_like_header(data_lines[0], delim):
            toks = data_lines[0].split(delim) if delim else data_lines[0].split()
            header_names = [t.strip() for t in toks][1:]
            data_lines = data_lines[1:]
        if decimal == "," and delim != ",":
            data_lines = [ln.replace(",", ".") for ln in data_lines]
        dat = np.genfromtxt(data_lines, delimiter=delim, dtype=float,
                            comments="#", invalid_raise=False)
        if dat.ndim == 1:
            dat = dat.reshape(-1, 1)
        dat = dat[~np.isnan(dat[:, 0])]
        if dat.size == 0 or dat.shape[1] < 2:
            raise ValueError("No numeric (X, Y) data found in the file.")
        arr = dat
        names = header_names if (header_names and len(header_names) >= dat.shape[1] - 1) else None

    x = arr[:, 0]
    spectra = []
    n_y = arr.shape[1] - 1
    for j in range(n_y):
        y = arr[:, j + 1]
        if n_y == 1:
            # one spectrum per file -> the filename is the clearest label
            name = base
        elif names:
            name = "{} · {}".format(base, names[j])
        else:
            name = "{} [{}]".format(base, j + 1)
        spectra.append((name, y))
    return x, spectra


class _Step(QtWidgets.QGroupBox):
    """A checkable, numbered pipeline step with a one-line hint."""

    changed = QtCore.pyqtSignal()

    def __init__(self, number, title, hint, parent=None):
        super().__init__("①②③④⑤⑥"[number - 1] + "  " + title, parent)
        self.setCheckable(True)
        self.setChecked(False)
        self.setToolTip(hint)
        self._form = QtWidgets.QFormLayout(self)
        self._form.setContentsMargins(7, 7, 7, 5)   # top room so the hint clears the title
        self._form.setSpacing(3)
        self._form.setLabelAlignment(QtCore.Qt.AlignLeft)
        self._form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        # wrap the field below its label when the row is too wide -> no horizontal
        # scrollbar / clipped controls in a narrow panel.
        self._form.setRowWrapPolicy(QtWidgets.QFormLayout.WrapLongRows)
        hint_lbl = QtWidgets.QLabel(hint)
        hint_lbl.setWordWrap(True)
        hint_lbl.setStyleSheet("color: #5a6472; font-size: 7pt;")
        self._form.addRow(hint_lbl)
        self.toggled.connect(lambda *_: self.changed.emit())

    def add_row(self, label, widget):
        self._form.addRow(label, widget)


class ResultPreviewDialog(QtWidgets.QDialog):
    """Clean, full-scale preview of the final (baseline-subtracted) spectrum,
    shown on its own so the user can review it before exporting."""

    def __init__(self, x, y, peaks, name, baseline_done, normalized, parent=None, fitted=None):
        super().__init__(parent)
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.peaks = peaks
        self.fitted = fitted or []
        self.name = name
        kind = "Baseline-subtracted spectrum" if baseline_done else "Processed spectrum"
        if normalized:
            kind += " (normalized)"
        self.setWindowTitle("Result preview — {}".format(name))
        self.resize(860, 580)

        layout = QtWidgets.QVBoxLayout(self)
        self.figure = Figure(tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        ax = self.figure.add_subplot(111)
        ax.plot(self.x, self.y, color="#1f6fb2", lw=1.4, label="spectrum")
        ax.axhline(0, color="#999999", lw=0.6)
        if self.fitted:
            model = np.zeros_like(self.x)
            for f in self.fitted:
                model = model + _lorentzian(self.x, f["pos"], f["h"], f["fwhm"])
            ax.plot(self.x, model, color="#2ca02c", lw=1.1, alpha=0.9, label="Lorentzian fit")
            ax.legend(loc="best", fontsize=8)
        if peaks:
            px = [p[0] for p in peaks]
            py = [p[1] for p in peaks]
            ax.plot(px, py, "v", color="#cc3344", ms=7)
            for pos, h in peaks:
                ax.annotate("{:.0f}".format(pos), (pos, h), textcoords="offset points",
                            xytext=(0, 6), ha="center", fontsize=7, color="#cc3344")
        ax.set_title("{}   ·   {} peaks".format(kind, len(peaks)))
        ax.set_xlabel("Raman shift / cm$^{-1}$")
        ax.set_ylabel("Intensity" + (" (normalized)" if normalized else " / a.u."))
        ax.grid(True, alpha=0.25)
        # tight y-margins so the result fills the view (the whole point of preview)
        ax.margins(x=0.01)

        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas)

        btns = QtWidgets.QHBoxLayout()
        b_data = QtWidgets.QPushButton("Export Data…")
        b_data.clicked.connect(self._export_data)
        b_fig = QtWidgets.QPushButton("Export Figure…")
        b_fig.clicked.connect(self._export_figure)
        b_close = QtWidgets.QPushButton("Close")
        b_close.clicked.connect(self.accept)
        btns.addWidget(b_data)
        btns.addWidget(b_fig)
        btns.addStretch(1)
        btns.addWidget(b_close)
        layout.addLayout(btns)

    def _export_data(self):
        fn, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export processed spectrum", "", "Text/CSV (*.csv *.txt)")
        if not fn:
            return
        if not os.path.splitext(fn)[1]:
            fn += ".csv"
        np.savetxt(fn, np.column_stack([self.x, self.y]), delimiter=",",
                   header="Raman shift (cm-1),Intensity", comments="")

    def _export_figure(self):
        fn, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export figure", "", "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if not fn:
            return
        self.figure.savefig(fn, dpi=200, bbox_inches="tight")


class SpectrumProcessor(QtWidgets.QMainWindow):
    """Guided, live-preview processing view for one or more spectra sharing an x axis."""

    def __init__(self, parent=None):
        super().__init__(parent)
        # library of all loaded spectra (each independent: its own x axis)
        self.library = []          # list of {name, x, y, source}
        # active spectrum (the one being processed)
        self.x_full = np.linspace(0, 1, 2)
        self.y_full = np.zeros_like(self.x_full)

        self.peaks = []            # list of (position, height)
        self.baseline_curve = None  # (x, baseline) to overlay, or None
        self.pre_baseline_y = None
        self.overlays = []         # [(name, x, y)] faint comparison traces
        self.xp = self.x_full
        self.yp = self.y_full

        # manual baseline anchors (step 4 "Edit anchors"): the baseline becomes a
        # monotone-cubic (PCHIP) interpolation through draggable anchor points.
        self.anchor_mode = False
        self.anchors_x = None
        self.anchors_y = None
        self._drag_idx = None
        self._hover_idx = None
        self._preserve_view = False
        self.fitted = []          # list of dicts {pos,h,fwhm,area} from Fit Peaks

        self.setWindowTitle("PyRaman Pro — Spectrum Processor")
        self._build_ui()
        self._update_empty_state()

    # -- spectra library ----------------------------------------------------
    def add_spectrum(self, name, x, y, source=""):
        """Add one spectrum to the library and the left list."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        n = min(len(x), len(y))
        self.library.append({"name": str(name), "x": x[:n], "y": y[:n], "source": source})
        item = QtWidgets.QListWidgetItem(str(name))
        item.setToolTip(source or str(name))
        self.lst_spectra.addItem(item)

    def load_files(self, paths):
        """Load one or more data files into the library; select the first new one."""
        first_new = len(self.library)
        errors = []
        for p in paths:
            try:
                x, spectra = load_spectra_file(p)
                for name, y in spectra:
                    self.add_spectrum(name, x, y, source=p)
            except Exception as exc:  # noqa: BLE001
                errors.append("{}: {}".format(os.path.basename(p), exc))
        if errors:
            QtWidgets.QMessageBox.warning(
                self, "Some files could not be loaded", "\n".join(errors))
        if len(self.library) > first_new:
            self.lst_spectra.setCurrentRow(first_new)  # triggers selection
        self._update_empty_state()

    # -- UI -----------------------------------------------------------------
    def _build_ui(self):
        self._build_menubar()

        central = QtWidgets.QWidget()
        outer = QtWidgets.QVBoxLayout(central)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # top action bar
        bar = QtWidgets.QToolBar()
        bar.setMovable(False)
        bar.setIconSize(QtCore.QSize(18, 18))
        act_auto = bar.addAction("Auto Analyze")
        act_auto.setToolTip("Run automatic baseline correction + peak detection, then fine-tune.")
        act_auto.triggered.connect(self.auto_analyze)
        if not _AUTORAMAN:
            act_auto.setEnabled(False)
        act_reset = bar.addAction("Reset")
        act_reset.triggered.connect(self.reset_steps)
        bar.addSeparator()
        act_preview = bar.addAction("Preview Result")
        act_preview.setToolTip("Open a clean, full-scale view of the baseline-subtracted spectrum before exporting.")
        act_preview.triggered.connect(self.preview_result)
        bar.addSeparator()
        act_exp_data = bar.addAction("Export Data")
        act_exp_data.triggered.connect(self.export_data)
        act_exp_fig = bar.addAction("Export Figure")
        act_exp_fig.triggered.connect(self.export_figure)
        bar.addSeparator()
        act_resetview = bar.addAction("Reset view")
        act_resetview.setToolTip("Zoom out to the full spectrum (scroll = zoom · drag = pan).")
        act_resetview.triggered.connect(self.reset_view)
        outer.addWidget(bar)

        # main split: spectra list | plot | processing panel
        split = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        split.addWidget(self._build_library_panel())

        plot_box = QtWidgets.QWidget()
        plot_layout = QtWidgets.QVBoxLayout(plot_box)
        plot_layout.setContentsMargins(0, 0, 0, 0)
        self.figure = Figure(figsize=(7, 5), tight_layout=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot(111)
        # interaction: mouse-wheel zoom (cursor-centred), left-drag pan, and
        # anchor dragging when Edit anchors is on. No matplotlib toolbar -- it
        # was bulky and most buttons went unused; scroll/drag is more intuitive
        # for zooming into anchors.
        self.canvas.mpl_connect("button_press_event", self._on_press)
        self.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.canvas.mpl_connect("button_release_event", self._on_release)
        self.canvas.mpl_connect("scroll_event", self._on_scroll)
        self._pan = None  # (x0px, y0px, xlim, ylim) while panning
        plot_layout.addWidget(self.canvas)
        hint = QtWidgets.QLabel("scroll = zoom · drag = pan · Edit anchors: drag = move, "
                                "double-click = add, right-click = delete")
        hint.setStyleSheet("color: #8a93a0; font-size: 8pt; padding: 2px 6px;")
        hint.setAlignment(QtCore.Qt.AlignHCenter)
        plot_layout.addWidget(hint)
        split.addWidget(plot_box)

        split.addWidget(self._build_panel())
        split.setStretchFactor(0, 0)
        split.setStretchFactor(1, 1)   # plot takes all extra width
        split.setStretchFactor(2, 0)
        split.setSizes([140, 1200, 230])
        self._split = split
        outer.addWidget(split)

        self.setCentralWidget(central)
        self.status = self.statusBar()
        self.status.showMessage("Open a spectrum file to begin (File ▸ Open files…).")

    def _build_menubar(self):
        mbar = self.menuBar()
        m_file = mbar.addMenu("File")
        act_open = m_file.addAction("Open files…", self.open_files)
        act_open.setShortcut("Ctrl+O")
        m_file.addSeparator()
        m_file.addAction("Export Data…", self.export_data)
        m_file.addAction("Export Figure…", self.export_figure)
        m_file.addSeparator()
        act_exit = m_file.addAction("Exit", self.close)
        act_exit.setShortcut("Ctrl+Q")

        m_tools = mbar.addMenu("Tools")
        m_tools.addAction("Fit peaks (FWHM / area)", self.fit_peaks)

        m_help = mbar.addMenu("Help")
        m_help.addAction("About", self._about)

    def _build_library_panel(self):
        box = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(box)
        v.setContentsMargins(8, 8, 6, 8)
        v.setSpacing(6)
        title = QtWidgets.QLabel("SPECTRA")
        title.setStyleSheet("font-weight: bold; letter-spacing: 1px; color: #33404f;")
        v.addWidget(title)
        self.lst_spectra = QtWidgets.QListWidget()
        self.lst_spectra.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.lst_spectra.currentRowChanged.connect(self._on_library_current_changed)
        self.lst_spectra.itemSelectionChanged.connect(self._on_library_selection_changed)
        v.addWidget(self.lst_spectra, 1)
        btn_open = QtWidgets.QPushButton("➕  Open files…")
        btn_open.clicked.connect(self.open_files)
        v.addWidget(btn_open)
        box.setMinimumWidth(160)
        return box

    def _update_empty_state(self):
        """Enable/disable processing when there is (no) data loaded."""
        has = len(self.library) > 0
        self.centralWidget().setEnabled(True)
        # the processing panel is meaningless with no data
        if hasattr(self, "panel_scroll"):
            self.panel_scroll.setEnabled(has)
        if not has:
            self.ax.clear()
            self.ax.text(0.5, 0.5, "Open a spectrum file to begin\n(File ▸ Open files…)",
                         ha="center", va="center", color="#8a93a0",
                         transform=self.ax.transAxes, fontsize=12)
            self.ax.set_xticks([]); self.ax.set_yticks([])
            self.canvas.draw_idle()

    def _build_panel(self):
        panel = QtWidgets.QScrollArea()
        panel.setWidgetResizable(True)
        panel.setFrameShape(QtWidgets.QFrame.NoFrame)
        panel.setMinimumWidth(200)
        panel.setMaximumWidth(280)
        # horizontal scrollbar should never be needed -- rows wrap instead
        panel.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.panel_scroll = panel
        inner = QtWidgets.QWidget()
        # NOTE: the global app stylesheet has `* { font-size: 10pt }`, and a Qt
        # stylesheet rule overrides widget.setFont(). So the panel font MUST be
        # set via a stylesheet here (which is more specific) to take effect.
        inner.setStyleSheet(
            "QWidget { font-size: 7pt; }"
            "QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit { min-height: 15px; }"
            "QPushButton { padding: 2px 6px; min-height: 15px; }"
            "QGroupBox { margin-top: 12px; }"           # room for the title
            "QGroupBox::title { font-size: 8pt; }")
        v = QtWidgets.QVBoxLayout(inner)
        v.setContentsMargins(6, 5, 6, 5)
        v.setSpacing(4)

        title = QtWidgets.QLabel("PROCESSING")
        title.setStyleSheet("font-weight: bold; letter-spacing: 1px; color: #33404f; font-size: 8pt;")
        v.addWidget(title)

        xr0, xr1 = float(self.x_full.min()), float(self.x_full.max())

        # 1. Range
        self.step_range = _Step(1, "Range (crop)", "Limit the spectrum to a wavenumber window.")
        self.sb_rmin = QtWidgets.QDoubleSpinBox(); self.sb_rmin.setRange(xr0, xr1); self.sb_rmin.setValue(xr0)
        self.sb_rmax = QtWidgets.QDoubleSpinBox(); self.sb_rmax.setRange(xr0, xr1); self.sb_rmax.setValue(xr1)
        for sb in (self.sb_rmin, self.sb_rmax):
            sb.setDecimals(0); sb.setSuffix(" cm⁻¹"); sb.valueChanged.connect(self._recompute)
        self.step_range.add_row("From", self.sb_rmin)
        self.step_range.add_row("To", self.sb_rmax)
        self.step_range.changed.connect(self._recompute)
        v.addWidget(self.step_range)

        # 2. Despike
        self.step_despike = _Step(2, "Despike", "Remove single-point cosmic-ray spikes.")
        self.sb_despike = QtWidgets.QDoubleSpinBox(); self.sb_despike.setRange(2, 15); self.sb_despike.setValue(5); self.sb_despike.setSingleStep(0.5)
        self.sb_despike.valueChanged.connect(self._recompute)
        self.step_despike.add_row("Threshold (σ)", self.sb_despike)
        self.step_despike.changed.connect(self._recompute)
        v.addWidget(self.step_despike)

        # 3. Smooth
        self.step_smooth = _Step(3, "Smooth", "Savitzky-Golay smoothing to reduce noise.")
        self.sb_win = QtWidgets.QSpinBox(); self.sb_win.setRange(3, 51); self.sb_win.setSingleStep(2); self.sb_win.setValue(7)
        self.sb_poly = QtWidgets.QSpinBox(); self.sb_poly.setRange(1, 5); self.sb_poly.setValue(3)
        self.sb_win.valueChanged.connect(self._recompute); self.sb_poly.valueChanged.connect(self._recompute)
        self.step_smooth.add_row("Window", self.sb_win)
        self.step_smooth.add_row("Poly order", self.sb_poly)
        self.step_smooth.changed.connect(self._recompute)
        v.addWidget(self.step_smooth)

        # 4. Baseline
        self.step_baseline = _Step(4, "Baseline", "Subtract the fluorescence background.")
        self.cmb_baseline = QtWidgets.QComboBox()
        for name in BASELINE_METHODS:
            self.cmb_baseline.addItem(name)
        self.cmb_baseline.currentIndexChanged.connect(self._baseline_method_changed)
        self.sl_strength = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.lbl_strength = QtWidgets.QLabel("Strength")
        self.sl_strength.valueChanged.connect(self._on_strength_changed)
        self.step_baseline.add_row("Method", self.cmb_baseline)
        self.step_baseline.add_row(self.lbl_strength, self.sl_strength)
        # manual anchor editing: turn the auto baseline into draggable points
        self.btn_anchors = QtWidgets.QPushButton("✎  Edit anchors")
        self.btn_anchors.setCheckable(True)
        self.btn_anchors.setToolTip(
            "Convert the baseline into draggable anchor points.\n"
            "Drag a point to reshape the baseline · double-click to add · "
            "right-click to delete.")
        self.btn_anchors.toggled.connect(self._toggle_anchor_mode)
        self.btn_reset_anchors = QtWidgets.QPushButton("Reset anchors")
        self.btn_reset_anchors.setToolTip("Re-seed the anchors from the automatic baseline.")
        self.btn_reset_anchors.clicked.connect(self._reset_anchors)
        self.btn_reset_anchors.setEnabled(False)
        self.btn_find = QtWidgets.QPushButton("Find anchors")
        self.btn_find.setToolTip(
            "Re-find anchors automatically on the low-curvature (background) "
            "points, like OriginLab's 2nd-derivative method.")
        self.btn_find.clicked.connect(self._find_anchors)
        self.btn_find.setEnabled(False)
        anchor_row = QtWidgets.QHBoxLayout()
        anchor_row.addWidget(self.btn_anchors)
        anchor_row.addWidget(self.btn_reset_anchors)
        anchor_row.addWidget(self.btn_find)
        self.step_baseline.add_row("Manual", anchor_row)
        # number of anchor handles. Dense by default (like OriginLab's ~60) so the
        # B-spline follows the broad background closely while skipping sharp peaks.
        self.sb_anchor_n = QtWidgets.QSpinBox()
        self.sb_anchor_n.setRange(4, 200)
        self.sb_anchor_n.setValue(40)
        self.sb_anchor_n.setToolTip("How many anchors to place. Changing this re-seeds them.")
        self.sb_anchor_n.valueChanged.connect(self._on_anchor_count_changed)
        self.step_baseline.add_row("Anchor count", self.sb_anchor_n)
        # interpolation through the anchors (OriginLab uses B-spline)
        self.cmb_interp = QtWidgets.QComboBox()
        self.cmb_interp.addItems(["B-spline", "PCHIP", "Linear"])
        self.cmb_interp.setToolTip("How the baseline is drawn through the anchors.")
        self.cmb_interp.currentIndexChanged.connect(self._recompute)
        self.step_baseline.add_row("Interpolation", self.cmb_interp)
        # flatness threshold: keep anchors only on points whose curvature is below
        # this percentile (lower = stricter = only the flattest background points).
        self.sb_anchor_thr = QtWidgets.QSpinBox()
        self.sb_anchor_thr.setRange(20, 95)
        self.sb_anchor_thr.setValue(70)
        self.sb_anchor_thr.setSuffix(" %")
        self.sb_anchor_thr.setToolTip(
            "Curvature threshold for auto-finding anchors: keep points flatter than "
            "this percentile (the background), excluding the curved peak regions.")
        self.sb_anchor_thr.valueChanged.connect(self._on_anchor_count_changed)
        self.step_baseline.add_row("Flatness", self.sb_anchor_thr)
        # vertical-only drag: dragging never shifts an anchor's wavenumber, which
        # makes manual baseline shaping much easier (x set by add/delete).
        self.cb_lockx = QtWidgets.QCheckBox("Lock anchor x (vertical drag)")
        self.cb_lockx.setChecked(True)
        self.cb_lockx.setToolTip("Dragging moves an anchor up/down only; its position on the x axis stays put.")
        self.step_baseline.add_row("", self.cb_lockx)
        # snap dropped anchors onto the true valley floor (keeps baseline under
        # the signal so the corrected spectrum doesn't dip negative)
        self.cb_snap = QtWidgets.QCheckBox("Snap anchors to valley")
        self.cb_snap.setChecked(True)
        self.cb_snap.setToolTip("When you ADD an anchor (double-click), lower it to the local "
                                "minimum nearby. Manual drags are kept exactly where you drop them.")
        self.step_baseline.add_row("", self.cb_snap)
        # clamp the corrected signal at 0 (the standard reference rests at zero)
        self.cb_clamp = QtWidgets.QCheckBox("Clamp result ≥ 0")
        self.cb_clamp.setChecked(False)
        self.cb_clamp.setToolTip("Set any negative corrected values to zero, like a clean reference spectrum.")
        self.cb_clamp.toggled.connect(self._recompute)
        self.step_baseline.add_row("", self.cb_clamp)
        # IS-Score: objective baseline-quality readout (Innocente et al., ACS
        # Omega 2026). Updates live whenever a baseline is applied.
        self.lbl_isscore = QtWidgets.QLabel("—")
        self.lbl_isscore.setToolTip(
            "IS-Score (Integrity Spectrum Score), 0–1, higher = better.\n"
            "Ground-truth-free baseline-quality metric (Innocente et al.,\n"
            "ACS Omega 2026). Green ≥ 0.8, amber ≥ 0.6, red below.")
        if _AUTORAMAN:
            self.step_baseline.add_row("IS-Score", self.lbl_isscore)
        self.btn_optimize = QtWidgets.QPushButton("Optimize (IS-Score)")
        self.btn_optimize.setToolTip(
            "Search this method's parameters to maximise the IS-Score, then set "
            "the controls to the best result.")
        self.btn_optimize.clicked.connect(self.optimize_baseline)
        if _AUTORAMAN:
            self.step_baseline.add_row("", self.btn_optimize)
        self.step_baseline.changed.connect(self._recompute)
        if not BASELINE_METHODS:
            self.step_baseline.setEnabled(False)
        v.addWidget(self.step_baseline)

        # 5. Normalize
        self.step_norm = _Step(5, "Normalize", "Scale intensity for comparison between spectra.")
        self.cmb_norm = QtWidgets.QComboBox(); self.cmb_norm.addItems(["Min-Max (0..1)", "Max = 1", "Area = 1"])
        self.cmb_norm.currentIndexChanged.connect(self._recompute)
        self.step_norm.add_row("Mode", self.cmb_norm)
        self.step_norm.changed.connect(self._recompute)
        v.addWidget(self.step_norm)

        # 6. Peaks
        self.step_peaks = _Step(6, "Peaks", "Detect bands. Higher % = fewer, stronger peaks.")
        self.sb_peak_pct = QtWidgets.QDoubleSpinBox(); self.sb_peak_pct.setRange(1, 40)
        self.sb_peak_pct.setValue(8); self.sb_peak_pct.setSingleStep(1); self.sb_peak_pct.setSuffix(" %")
        self.sb_peak_pct.valueChanged.connect(self._recompute)
        self.sb_peak_pct.setToolTip("Minimum peak height/prominence as a percent of the strongest band.")
        self.step_peaks.add_row("Min height/prom %", self.sb_peak_pct)
        self.btn_fit = QtWidgets.QPushButton("Fit peaks  (FWHM / area)")
        self.btn_fit.setToolTip("Fit each detected band with a Lorentzian to get FWHM and area "
                                "(overlapping bands fitted jointly).")
        self.btn_fit.clicked.connect(self.fit_peaks)
        self.step_peaks.add_row("", self.btn_fit)
        self.step_peaks.changed.connect(self._recompute)
        if not _AUTORAMAN:
            self.step_peaks.setEnabled(False)
        v.addWidget(self.step_peaks)

        # peak table
        self.tbl_peaks = QtWidgets.QTableWidget(0, 4)
        self.tbl_peaks.setHorizontalHeaderLabels(["Position cm⁻¹", "Height", "FWHM", "Area"])
        self.tbl_peaks.horizontalHeader().setStretchLastSection(True)
        self.tbl_peaks.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.tbl_peaks.setMaximumHeight(220)
        v.addWidget(self.tbl_peaks)

        v.addStretch(1)
        self._baseline_method_changed()
        panel.setWidget(inner)
        return panel

    # -- behaviour ----------------------------------------------------------
    def _baseline_method_changed(self):
        if not BASELINE_METHODS:
            return
        spec = BASELINE_METHODS[self.cmb_baseline.currentText()]
        lo, hi, default, label, _ = KNOB_SPEC[spec["knob"]]
        self.sl_strength.blockSignals(True)
        self.sl_strength.setRange(lo, hi)
        self.sl_strength.setValue(default)
        self.sl_strength.blockSignals(False)
        self._update_strength_label()
        if self.step_baseline.isChecked():
            self._recompute()

    def _update_strength_label(self):
        """Show the live parameter value next to the strength slider."""
        if not BASELINE_METHODS:
            return
        spec = BASELINE_METHODS[self.cmb_baseline.currentText()]
        knob = spec["knob"]
        val = KNOB_SPEC[knob][4](self.sl_strength.value())
        if knob == "lam":
            self.lbl_strength.setText("Strength  λ={:.0e}".format(val))
        elif knob == "poly_order":
            self.lbl_strength.setText("Poly order  {}".format(val))
        else:
            self.lbl_strength.setText("Ball radius  {}".format(val))

    def _on_strength_changed(self):
        self._update_strength_label()
        self._recompute()

    def open_files(self):
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(
            self, "Open spectrum file(s)", "",
            "Spectra (*.csv *.txt *.dat *.asc *.xlsx *.xls);;All files (*)")
        if paths:
            self.load_files(paths)

    def _on_library_current_changed(self, row):
        """Make the current row the active (processed) spectrum."""
        if 0 <= row < len(self.library):
            entry = self.library[row]
            self.x_full = entry["x"]
            self.y_full = entry["y"]
            self.active_name = entry["name"]
            # leaving anchor mode avoids stale anchors from another spectrum
            if self.btn_anchors.isChecked():
                self.btn_anchors.blockSignals(True)
                self.btn_anchors.setChecked(False)
                self.btn_anchors.blockSignals(False)
                self.anchor_mode = False
            self.setWindowTitle("PyRaman Pro — " + entry["name"])
            self._recompute()

    def _on_library_selection_changed(self):
        """Collect any extra-selected spectra as faint comparison overlays."""
        rows = sorted(i.row() for i in self.lst_spectra.selectedIndexes())
        cur = self.lst_spectra.currentRow()
        self.overlays = [
            (self.library[r]["name"], self.library[r]["x"], self.library[r]["y"])
            for r in rows if r != cur and 0 <= r < len(self.library)
        ]
        self._redraw()

    def optimize_baseline(self):
        """Tune the current method's parameters to maximise the IS-Score."""
        if not _AUTORAMAN or not self.library:
            return
        x, y = self._preprocessed()
        spec = BASELINE_METHODS[self.cmb_baseline.currentText()]
        try:
            res = ParameterOptimizer().optimize(spec["fn"], x, y, self._baseline_param()[1])
        except Exception as exc:  # noqa: BLE001
            self.status.showMessage("Optimize failed: {}".format(exc), 6000)
            return
        y_corr, baseline, params, score = res
        if y_corr is None:
            self.status.showMessage("Optimize: no valid baseline found.", 5000)
            return
        # reflect the best knob value back into the strength slider
        knob = spec["knob"]
        if knob in params:
            val = params[knob]
            if knob == "lam":
                slider_val = int(round(10.0 * np.log10(val)))
            else:
                slider_val = int(val)
            self.sl_strength.blockSignals(True)
            self.sl_strength.setValue(slider_val)
            self.sl_strength.blockSignals(False)
        if not self.step_baseline.isChecked():
            self.step_baseline.setChecked(True)
        self._recompute()
        self.status.showMessage("Optimized {}: IS-Score {:.3f}".format(
            self.cmb_baseline.currentText(), score), 8000)

    def _update_isscore(self, x, corrected):
        """Update the IS-Score readout for the current baseline correction."""
        if not _AUTORAMAN or is_score is None:
            return
        try:
            s = is_score(self.pre_baseline_y, corrected, x)
        except Exception:  # noqa: BLE001
            self.lbl_isscore.setText("—"); return
        color = "#2e8b57" if s >= 0.8 else ("#c77b1e" if s >= 0.6 else "#c0392b")
        self.lbl_isscore.setText("{:.3f}".format(s))
        self.lbl_isscore.setStyleSheet("font-weight: bold; color: {};".format(color))

    def _about(self):
        QtWidgets.QMessageBox.about(
            self, "About PyRaman Pro",
            "PyRaman Pro — Spectrum Processor\n\n"
            "Automatic baseline correction is scored with the IS-Score\n"
            "(Integrity Spectrum Score), a re-implementation of:\n\n"
            "  Innocente et al., \"Automated Baseline Correction Evaluation\n"
            "  Score for Raman Spectroscopy\", ACS Omega 2026, 11, 25057.\n"
            "  DOI: 10.1021/acsomega.5c09870")

    def _baseline_param(self):
        spec = BASELINE_METHODS[self.cmb_baseline.currentText()]
        knob = spec["knob"]
        value = KNOB_SPEC[knob][4](self.sl_strength.value())
        return spec["fn"], {knob: value}

    def _preprocessed(self):
        """Return (x, y) after the pre-baseline steps: range, despike, smooth."""
        x = self.x_full.copy()
        y = self.y_full.copy()
        if self.step_range.isChecked():
            lo, hi = sorted([self.sb_rmin.value(), self.sb_rmax.value()])
            m = (x >= lo) & (x <= hi)
            if m.sum() >= 5:
                x, y = x[m], y[m]
        if self.step_despike.isChecked():
            y = despike(y, self.sb_despike.value())
        if self.step_smooth.isChecked():
            w = self.sb_win.value()
            if w % 2 == 0:
                w += 1
            po = self.sb_poly.value()
            if w > po and w <= len(y):
                y = savgol_filter(y, w, po)
        return x, y

    def _anchor_baseline(self, x):
        """Baseline drawn through the current anchors with the selected method
        (B-spline like OriginLab, PCHIP, or Linear)."""
        order = np.argsort(self.anchors_x)
        ax_, ay_ = self.anchors_x[order], self.anchors_y[order]
        # de-duplicate identical x (interpolators need strictly increasing x)
        keep = np.concatenate(([True], np.diff(ax_) > 0))
        ax_, ay_ = ax_[keep], ay_[keep]
        if len(ax_) < 2:
            return np.full_like(x, ay_[0] if len(ay_) else 0.0)
        method = self.cmb_interp.currentText() if hasattr(self, "cmb_interp") else "B-spline"
        if method == "Linear":
            return np.interp(x, ax_, ay_)
        if method == "B-spline" and len(ax_) >= 4:
            try:
                from scipy.interpolate import make_interp_spline
                return make_interp_spline(ax_, ay_, k=3)(x)
            except Exception:  # noqa: BLE001 - fall back to PCHIP
                pass
        from scipy.interpolate import PchipInterpolator
        return PchipInterpolator(ax_, ay_, extrapolate=True)(x)

    # -- manual baseline anchors --------------------------------------------
    def _seed_anchors(self, n=None):
        """Auto-find anchors on the LOW-CURVATURE (background) points, the way
        OriginLab's "1st & 2nd Derivative (zeroes)" method does: baseline regions
        have small curvature, peaks have large curvature. We SG-smooth, take the
        2nd derivative, keep the points whose |curvature| is below the Flatness
        percentile (so peaks are excluded), spread ``n`` anchors evenly across
        those, and snap each onto the local minimum so the baseline sits under
        the signal. Dense anchors + a B-spline then follow the broad background
        while skipping the sharp Raman peaks.
        """
        if n is None:
            n = self.sb_anchor_n.value()
        x, y = self._preprocessed()
        npts = len(x)
        if npts < 4:
            self.anchors_x = x.astype(float)
            self.anchors_y = y.astype(float)
            return
        # smoothed 2nd derivative -> curvature
        w = max(7, (npts // 50) | 1)
        w = min(w, npts if npts % 2 else npts - 1)
        ys = savgol_filter(y, w, 3) if npts > w else y.copy()
        curv = np.abs(np.gradient(np.gradient(ys)))
        thr_pct = self.sb_anchor_thr.value() if hasattr(self, "sb_anchor_thr") else 70
        thr = np.percentile(curv, thr_pct)
        cand = np.where(curv <= thr)[0]              # low-curvature = background
        if cand.size < 2:
            cand = np.arange(npts)
        cand = np.unique(np.concatenate(([0], cand, [npts - 1])))
        # spread n anchors evenly across x, each snapped to nearest background pt
        targets = np.linspace(0, npts - 1, max(2, n))
        idx = np.unique([int(cand[np.argmin(np.abs(cand - t))]) for t in targets])
        # lower each anchor onto the true local minimum nearby (under the signal)
        hw = max(2, int(0.008 * npts))
        ay = np.array([float(np.min(y[max(0, i - hw): i + hw + 1])) for i in idx])
        self.anchors_x = x[idx].astype(float)
        self.anchors_y = ay

    def _find_anchors(self):
        """Re-run automatic anchor finding at the current settings."""
        if self.anchor_mode:
            self._seed_anchors()
            self._recompute()
            self.status.showMessage(
                "Found {} anchors on the background (Flatness {}%).".format(
                    len(self.anchors_x), self.sb_anchor_thr.value()), 6000)

    def _on_anchor_count_changed(self):
        """Re-seed anchors at the new density (only meaningful while editing)."""
        if self.anchor_mode:
            self._seed_anchors()
            self._recompute()

    def _toggle_anchor_mode(self, on):
        self.anchor_mode = bool(on)
        self.btn_reset_anchors.setEnabled(self.anchor_mode)
        self.btn_find.setEnabled(self.anchor_mode)
        # method/strength don't apply while editing anchors by hand
        self.cmb_baseline.setEnabled(not self.anchor_mode)
        self.sl_strength.setEnabled(not self.anchor_mode)
        if self.anchor_mode:
            if not self.step_baseline.isChecked():
                self.step_baseline.setChecked(True)  # triggers a recompute
            self._seed_anchors()
            self.status.showMessage(
                "Anchor editing on: anchors seeded on the valleys · drag up/down "
                "· double-click to add · right-click to delete.", 9000)
        self._recompute()

    def _reset_anchors(self):
        if self.anchor_mode:
            self._seed_anchors()
            self._recompute()

    def _nearest_anchor(self, event, tol_px=18):
        if self.anchors_x is None or len(self.anchors_x) == 0:
            return None
        pts = self.ax.transData.transform(np.column_stack([self.anchors_x, self.anchors_y]))
        d = np.hypot(pts[:, 0] - event.x, pts[:, 1] - event.y)
        i = int(np.argmin(d))
        return i if d[i] <= tol_px else None

    def _on_press(self, event):
        if event.inaxes != self.ax or event.xdata is None:
            return
        if self.anchor_mode:
            i = self._nearest_anchor(event)
            if event.button == 3:  # right-click -> delete nearest anchor
                if i is not None:
                    self._delete_anchor(i)
                return
            if i is not None:                     # grab an anchor to drag
                self._drag_idx = i
                return
            if event.dblclick:                    # double-click empty -> add
                self._add_anchor(event.xdata, event.ydata)
                return
        # otherwise a left-drag pans the view (works in any mode, empty space)
        if event.button == 1:
            self._pan = (event.x, event.y, self.ax.get_xlim(), self.ax.get_ylim())

    def _on_motion(self, event):
        # panning the view
        if self._pan is not None and event.x is not None:
            x0, y0, xlim, ylim = self._pan
            inv = self.ax.transData.inverted()
            p0 = inv.transform((x0, y0))
            p1 = inv.transform((event.x, event.y))
            dx, dy = p0[0] - p1[0], p0[1] - p1[1]
            self.ax.set_xlim(xlim[0] + dx, xlim[1] + dx)
            self.ax.set_ylim(ylim[0] + dy, ylim[1] + dy)
            self.canvas.draw_idle()
            return
        # not dragging an anchor: hover feedback so anchors are easy to grab
        if self._drag_idx is None:
            if not self.anchor_mode or event.inaxes != self.ax or event.xdata is None:
                return
            i = self._nearest_anchor(event)
            if i != self._hover_idx:
                self._hover_idx = i
                cur = QtCore.Qt.PointingHandCursor if i is not None else QtCore.Qt.ArrowCursor
                self.canvas.setCursor(QtGui.QCursor(cur))
            return
        if event.inaxes != self.ax or event.xdata is None:
            return
        # vertical-only drag (default) keeps the anchor's wavenumber fixed
        if not self.cb_lockx.isChecked():
            self.anchors_x[self._drag_idx] = float(event.xdata)
        self.anchors_y[self._drag_idx] = float(event.ydata)
        self._preserve_view = True
        self._recompute()
        self._preserve_view = False

    def _on_release(self, event):
        if self._pan is not None:
            self._pan = None
            return
        if self._drag_idx is not None:
            self._drag_idx = None
            # Keep the anchor exactly where the user dropped it -- do NOT snap it
            # back to the valley (that was pulling dragged anchors back down).
            # "Snap to valley" applies only to newly added anchors.
            self._recompute()  # final pass (peaks etc.)

    def _on_scroll(self, event):
        """Mouse-wheel zoom, centred on the cursor. Shift+wheel zooms Y."""
        if event.inaxes != self.ax or event.xdata is None:
            return
        base = 1.25
        scale = (1.0 / base) if event.button == "up" else base   # up = zoom in
        zoom_y = bool(event.key and "shift" in event.key)
        if zoom_y:
            cur, c = self.ax.get_ylim(), event.ydata
        else:
            cur, c = self.ax.get_xlim(), event.xdata
        width = cur[1] - cur[0]
        rel = (c - cur[0]) / width if width else 0.5
        new_w = width * scale
        lo, hi = c - new_w * rel, c + new_w * (1 - rel)
        if zoom_y:
            self.ax.set_ylim(lo, hi)
        else:
            self.ax.set_xlim(lo, hi)
        self.canvas.draw_idle()

    def reset_view(self):
        """Zoom back out to the full spectrum."""
        if self.library:
            self._redraw()   # replots and autoscales to the data

    def _snap_anchor(self, i, half_width_frac=0.02):
        """Lower an anchor onto the local minimum of the (pre-baseline) signal
        nearby, so the baseline rests under the spectrum rather than through it."""
        if not self.cb_snap.isChecked() or self.anchors_x is None:
            return
        x, y = self._preprocessed()
        xi = self.anchors_x[i]
        hw = half_width_frac * (float(x.max()) - float(x.min()))
        m = (x >= xi - hw) & (x <= xi + hw)
        if m.any():
            yw = y[m]
            self.anchors_y[i] = float(np.min(yw))

    def _add_anchor(self, x, y):
        self.anchors_x = np.append(self.anchors_x, float(x))
        self.anchors_y = np.append(self.anchors_y, float(y))
        self._snap_anchor(len(self.anchors_x) - 1)
        self._preserve_view = True
        self._recompute()
        self._preserve_view = False

    def _delete_anchor(self, i):
        if self.anchors_x is None or len(self.anchors_x) <= 2:
            self.status.showMessage("Keep at least 2 anchors.", 4000)
            return
        self.anchors_x = np.delete(self.anchors_x, i)
        self.anchors_y = np.delete(self.anchors_y, i)
        self._preserve_view = True
        self._recompute()
        self._preserve_view = False

    def _recompute(self):
        if not self.library:
            self._update_empty_state()
            return
        self.baseline_curve = None
        self.pre_baseline_y = None

        x, y = self._preprocessed()

        if self.step_baseline.isChecked() and BASELINE_METHODS:
            try:
                self.pre_baseline_y = y.copy()
                if self.anchor_mode and self.anchors_x is not None and len(self.anchors_x) >= 2:
                    bl = self._anchor_baseline(x)
                else:
                    fn, params = self._baseline_param()
                    _, bl = fn(x, y, **params)
                    bl = np.asarray(bl, dtype=float).flatten()
                self.baseline_curve = (x.copy(), bl)
                corrected = y - bl
                self._update_isscore(x, corrected)  # objective quality readout
                y = corrected
                if self.cb_clamp.isChecked():
                    y = np.clip(y, 0.0, None)  # rest at zero like a clean reference
            except Exception as exc:  # noqa: BLE001
                self.status.showMessage("Baseline failed: {}".format(exc), 6000)
        elif _AUTORAMAN:
            self.lbl_isscore.setText("—")
            self.lbl_isscore.setStyleSheet("")

        if self.step_norm.isChecked():
            y = normalize(y, self.cmb_norm.currentText())

        self.xp, self.yp = x, y
        self.fitted = []  # any change invalidates a prior peak fit

        self.peaks = []
        if self.step_peaks.isChecked() and _AUTORAMAN:
            # 3-sigma noise gate keeps flat/noisy spectra from exploding; the
            # relative gate then requires a band to be at least N% as strong as
            # the tallest -- by BOTH prominence and height -- which removes the
            # small bumps in band-free regions (a flat reference has none there).
            detected = AutoPeaks(prominence_sigma=3.0, min_confidence=0.0).detect(x, y)
            if detected:
                pct = self.sb_peak_pct.value() / 100.0
                max_prom = max(p.prominence for p in detected)
                max_h = max(p.height for p in detected)
                detected = [p for p in detected
                            if p.prominence >= pct * max_prom and p.height >= pct * max_h]
            self.peaks = [(p.position, p.height) for p in detected]

        self._redraw()
        self._fill_peak_table()

    def _redraw(self):
        if self._preserve_view:
            xlim, ylim = self.ax.get_xlim(), self.ax.get_ylim()
        self.ax.clear()
        normalized = self.step_norm.isChecked()
        # context: raw (or pre-baseline) trace
        if self.step_baseline.isChecked() and self.pre_baseline_y is not None and not normalized:
            self.ax.plot(self.xp, self.pre_baseline_y, color="#c7ced8", lw=1.0, label="before baseline")
            if self.baseline_curve is not None:
                self.ax.plot(self.baseline_curve[0], self.baseline_curve[1],
                             color="#e08b3c", lw=1.2, ls="--", label="baseline")
                # draggable anchor handles
                if self.anchor_mode and self.anchors_x is not None and not normalized:
                    self.ax.plot(self.anchors_x, self.anchors_y, "o", color="#e08b3c",
                                 ms=9, mec="#7a4a16", mew=1.2, zorder=6,
                                 label="anchors (drag)")
        elif not normalized:
            self.ax.plot(self.x_full, self.y_full, color="#d3d9e0", lw=1.0, label="raw")
        # comparison overlays (other selected spectra, raw, faint)
        for name, ox, oy in self.overlays:
            self.ax.plot(ox, oy, lw=0.9, alpha=0.55, label=name)
        # processed
        self.ax.plot(self.xp, self.yp, color="#1f6fb2", lw=1.4, label="processed")
        # fitted model overlay (sum of Lorentzians)
        if self.fitted:
            model = np.zeros_like(self.xp)
            for f in self.fitted:
                model = model + _lorentzian(self.xp, f["pos"], f["h"], f["fwhm"])
            self.ax.plot(self.xp, model, color="#2ca02c", lw=1.1, alpha=0.9,
                         label="fit ({})".format(len(self.fitted)))
        # peaks
        if self.peaks:
            px = [p[0] for p in self.peaks]
            py = [p[1] for p in self.peaks]
            self.ax.plot(px, py, "v", color="#cc3344", ms=7, label="peaks ({})".format(len(self.peaks)))
            for pos, h in self.peaks[:20]:
                self.ax.annotate("{:.0f}".format(pos), (pos, h), textcoords="offset points",
                                 xytext=(0, 6), ha="center", fontsize=7, color="#cc3344")
        self.ax.set_xlabel("Raman shift / cm$^{-1}$")
        self.ax.set_ylabel("Intensity" + (" (normalized)" if normalized else " / a.u."))
        self.ax.legend(loc="best", fontsize=8, framealpha=0.85)
        self.ax.grid(True, alpha=0.25)
        if self._preserve_view:
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        self.canvas.draw_idle()

    def _fill_peak_table(self):
        if self.fitted:
            rows = sorted(self.fitted, key=lambda f: f["pos"])
            self.tbl_peaks.setRowCount(len(rows))
            for r, f in enumerate(rows):
                vals = ["{:.1f}".format(f["pos"]), "{:.3g}".format(f["h"]),
                        "{:.1f}".format(f["fwhm"]), "{:.3g}".format(f["area"])]
                for c, val in enumerate(vals):
                    self.tbl_peaks.setItem(r, c, QtWidgets.QTableWidgetItem(val))
        else:
            self.tbl_peaks.setRowCount(len(self.peaks))
            for r, (pos, h) in enumerate(self.peaks):
                self.tbl_peaks.setItem(r, 0, QtWidgets.QTableWidgetItem("{:.1f}".format(pos)))
                self.tbl_peaks.setItem(r, 1, QtWidgets.QTableWidgetItem("{:.3g}".format(h)))
                self.tbl_peaks.setItem(r, 2, QtWidgets.QTableWidgetItem(""))
                self.tbl_peaks.setItem(r, 3, QtWidgets.QTableWidgetItem(""))

    # -- actions ------------------------------------------------------------
    def auto_analyze(self):
        if not _AUTORAMAN:
            return
        if not self.library:
            self.status.showMessage("Open a spectrum first (File ▸ Open files…).", 4000)
            return
        report = AutoAnalyzer().analyze(self.x_full, self.y_full)
        if not report.ok:
            self.status.showMessage("Auto: baseline failed -- adjust the Baseline step manually.", 8000)
            return
        # leave manual-anchor mode so the auto baseline is used
        if self.btn_anchors.isChecked():
            self.btn_anchors.setChecked(False)  # -> _toggle_anchor_mode(False)
        # reflect the auto choices in the controls, then recompute
        self.step_smooth.setChecked(True)
        # match the chosen method by underlying function identity (the auto panel
        # uses long names, our dropdown uses short labels, but both call the same
        # function objects).
        idx = -1
        try:
            from autoRaman.auto_baseline import DEFAULT_METHOD_PANEL
            auto_fn = DEFAULT_METHOD_PANEL.get(report.baseline.method, (None,))[0]
            for i, label in enumerate(BASELINE_METHODS):
                if BASELINE_METHODS[label]["fn"] is auto_fn:
                    idx = i
                    break
        except Exception:  # noqa: BLE001
            idx = self.cmb_baseline.findText(report.baseline.method)
        if idx >= 0:
            self.cmb_baseline.blockSignals(True)
            self.cmb_baseline.setCurrentIndex(idx)
            self.cmb_baseline.blockSignals(False)
            self._baseline_method_changed()
        # push the auto-selected lambda into the strength slider so the manual
        # control reflects (and can fine-tune from) what Auto chose
        lam = report.baseline.params.get("lam")
        if lam:
            self.sl_strength.blockSignals(True)
            self.sl_strength.setValue(int(round(10.0 * np.log10(lam))))
            self.sl_strength.blockSignals(False)
        self.step_baseline.setChecked(True)
        self.step_peaks.setChecked(True)
        self._recompute()
        self.status.showMessage(report.summary().splitlines()[0] +
                                "  |  {} peaks".format(len(self.peaks)), 10000)

    def fit_peaks(self):
        """Fit each detected band with a Lorentzian (overlapping bands jointly)
        to report position, height, FWHM and area -- like OriginLab's Peak
        Analyzer. Operates on the current corrected spectrum."""
        if not self.peaks or not self.step_peaks.isChecked():
            self.status.showMessage("Enable the Peaks step and detect peaks first.", 5000)
            return
        x = np.asarray(self.xp, dtype=float)
        y = np.asarray(self.yp, dtype=float)
        peaks = sorted(self.peaks, key=lambda p: p[0])

        # group bands whose neighbourhoods overlap so they are fitted together
        merge = 60.0  # cm^-1
        groups, cur = [], [peaks[0]]
        for p in peaks[1:]:
            if p[0] - cur[-1][0] <= merge:
                cur.append(p)
            else:
                groups.append(cur)
                cur = [p]
        groups.append(cur)

        fitted, failed = [], 0
        for g in groups:
            xs = [p[0] for p in g]
            lo, hi = min(xs) - 40.0, max(xs) + 40.0
            m = (x >= lo) & (x <= hi)
            if m.sum() < len(g) * 3 + 2:
                continue
            xw, yw = x[m], y[m]
            p0, blo, bhi = [], [], []
            for px, ph in g:
                p0 += [px, max(ph, 1e-6), 15.0]
                blo += [px - 25.0, 0.0, 2.0]
                bhi += [px + 25.0, np.inf, 200.0]
            p0 += [0.0, 0.0]; blo += [-np.inf, -np.inf]; bhi += [np.inf, np.inf]
            try:
                popt, _ = curve_fit(_multi_lorentzian, xw, yw, p0=p0,
                                    bounds=(blo, bhi), maxfev=8000)
            except Exception:  # noqa: BLE001
                failed += 1
                continue
            for i in range(0, len(g) * 3, 3):
                xc, h, fw = popt[i], popt[i + 1], popt[i + 2]
                fitted.append({"pos": float(xc), "h": float(h), "fwhm": float(fw),
                               "area": float(np.pi * h * fw / 2.0)})

        if not fitted:
            self.status.showMessage("Peak fitting did not converge.", 6000)
            return
        self.fitted = fitted
        self._redraw()
        self._fill_peak_table()
        msg = "Fitted {} peaks (FWHM/area in table).".format(len(fitted))
        if failed:
            msg += " {} group(s) failed to converge.".format(failed)
        self.status.showMessage(msg, 8000)

    def reset_steps(self):
        if self.btn_anchors.isChecked():
            self.btn_anchors.setChecked(False)
        for step in (self.step_range, self.step_despike, self.step_smooth,
                     self.step_baseline, self.step_norm, self.step_peaks):
            step.setChecked(False)
        self.status.showMessage("Reset to raw spectrum.", 4000)
        self._recompute()

    def preview_result(self):
        """Show the final subtracted spectrum on its own, full-scale, with an
        Export option right there."""
        if not self.library:
            self.status.showMessage("Open a spectrum first.", 4000)
            return
        dlg = ResultPreviewDialog(
            self.xp, self.yp, self.peaks,
            getattr(self, "active_name", "Spectrum"),
            baseline_done=self.step_baseline.isChecked() and self.baseline_curve is not None,
            normalized=self.step_norm.isChecked(),
            parent=self, fitted=self.fitted,
        )
        dlg.exec_()

    def export_data(self):
        fn, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export processed spectrum", "", "Text/CSV (*.csv *.txt)")
        if not fn:
            return
        if not os.path.splitext(fn)[1]:
            fn += ".csv"
        np.savetxt(fn, np.column_stack([self.xp, self.yp]), delimiter=",",
                   header="Raman shift (cm-1),Intensity", comments="")
        self.status.showMessage("Saved {}".format(fn), 5000)

    def export_figure(self):
        fn, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export figure", "", "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if not fn:
            return
        self.figure.savefig(fn, dpi=200, bbox_inches="tight")
        self.status.showMessage("Saved {}".format(fn), 5000)
