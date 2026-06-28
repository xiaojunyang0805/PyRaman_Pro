"""
PyRaman Pro -- launcher.

A simple, spectrum-centric Raman analysis app: it opens straight into the
Spectrum Processor (spectrumProcessor.SpectrumProcessor) -- a Spectra list on the
left, the plot in the centre, and a step-by-step Processing panel on the right
(with a live IS-Score for baseline quality).

This module only sets up the QApplication look-and-feel and launches that window;
all functionality lives in spectrumProcessor.py + the autoRaman package.
"""

import os
import sys

from PyQt5 import QtCore, QtGui, QtWidgets

from spectrumProcessor import SpectrumProcessor


# Clean, flat light theme. An explicit sans-serif family is important: without
# it, Qt on Chinese Windows falls back to a CJK font whose Latin glyphs render
# wide and serif, which clipped the panel text. "Microsoft YaHei UI" is listed
# after Segoe UI so CJK characters still render.
APP_STYLESHEET = """
* { font-family: "Segoe UI", "Microsoft YaHei UI", "Microsoft YaHei", sans-serif; font-size: 10pt; }
QMenuBar::item { padding: 4px 10px; }
QMenu::item { padding: 5px 24px; }
QTreeWidget::item, QListWidget::item { padding: 3px 2px; }
QTableWidget { gridline-color: #c8c8c8; }
QHeaderView::section {
    font-weight: bold;
    padding: 4px 6px;
    background-color: #eef2f7;
    border: 1px solid #c8c8c8;
}
QPushButton { padding: 5px 14px; min-height: 22px; }
QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit { min-height: 24px; padding: 2px 4px; }
QTabBar::tab { padding: 6px 14px; }
QToolTip { font-size: 9pt; }
QStatusBar { font-size: 9pt; }
QToolBar { background: #e9edf2; border: none; padding: 4px; spacing: 6px; }
QToolBar QToolButton { padding: 5px 12px; border-radius: 4px; }
QToolBar QToolButton:hover { background: #d6deea; }
QGroupBox {
    font-weight: bold;
    border: 1px solid #d2d9e2;
    border-radius: 6px;
    margin-top: 10px;
    background: #fbfcfe;
}
QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 5px; color: #33404f; }
QGroupBox::indicator { width: 16px; height: 16px; }
"""


def apply_app_style(app):
    """Fusion style + flat light palette + an explicit sans-serif UI font.

    High-DPI scaling is enabled in main() before the QApplication is created.
    """
    app.setStyle("Fusion")

    palette = QtGui.QPalette()
    palette.setColor(QtGui.QPalette.Window, QtGui.QColor("#f4f6f9"))
    palette.setColor(QtGui.QPalette.WindowText, QtGui.QColor("#222831"))
    palette.setColor(QtGui.QPalette.Base, QtGui.QColor("#ffffff"))
    palette.setColor(QtGui.QPalette.AlternateBase, QtGui.QColor("#eef2f7"))
    palette.setColor(QtGui.QPalette.Text, QtGui.QColor("#222831"))
    palette.setColor(QtGui.QPalette.Button, QtGui.QColor("#e9edf2"))
    palette.setColor(QtGui.QPalette.ButtonText, QtGui.QColor("#222831"))
    palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor("#2f7fd1"))
    palette.setColor(QtGui.QPalette.HighlightedText, QtGui.QColor("#ffffff"))
    palette.setColor(QtGui.QPalette.ToolTipBase, QtGui.QColor("#ffffff"))
    palette.setColor(QtGui.QPalette.ToolTipText, QtGui.QColor("#222831"))
    app.setPalette(palette)

    # pick the first available clean sans-serif family (Latin), CJK falls back
    # via the stylesheet font-family list above.
    families = set(QtGui.QFontDatabase().families())
    for fam in ("Segoe UI", "Arial", "Helvetica"):
        if fam in families:
            app.setFont(QtGui.QFont(fam, 10))
            break
    else:
        base = app.font()
        base.setPointSize(10)
        app.setFont(base)
    app.setStyleSheet(APP_STYLESHEET)


def main():
    # High-DPI scaling must be set before the QApplication is instantiated.
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)
    app = QtWidgets.QApplication(sys.argv)
    apply_app_style(app)

    icon_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                             "Icons", "PyRAMAN_logo.png")
    if os.path.exists(icon_path):
        app.setWindowIcon(QtGui.QIcon(icon_path))

    win = SpectrumProcessor()
    if os.path.exists(icon_path):
        win.setWindowIcon(QtGui.QIcon(icon_path))
    win.showMaximized()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
