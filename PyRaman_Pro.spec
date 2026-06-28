# -*- mode: python ; coding: utf-8 -*-

import os
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

block_cipher = None

# Get the src directory
src_dir = os.path.join(os.getcwd(), 'src')

# Collect all data files
datas = []
# Add Icons folder
icons_path = os.path.join(src_dir, 'Icons')
if os.path.exists(icons_path):
    datas.append((icons_path, 'Icons'))

# Collect hidden imports
hiddenimports = [
    'numpy',
    'matplotlib',
    'scipy',
    'pandas',
    'openpyxl',
    'pybaselines',
    # AutoRaman is imported via a guarded `try: from autoRaman import ...`,
    # which PyInstaller's static analysis can miss -> list it explicitly.
    'autoRaman',
    'autoRaman.auto_analyzer',
    'autoRaman.auto_baseline',
    'autoRaman.auto_peaks',
    'autoRaman.quality_metrics',
    'autoRaman.is_score',
    'matplotlib.backends.backend_qt5agg',
    'PyQt5.QtCore',
    'PyQt5.QtGui',
    'PyQt5.QtWidgets',
]

# Collect submodules for packages that might need them
hiddenimports += collect_submodules('scipy')
hiddenimports += collect_submodules('matplotlib')

a = Analysis(
    ['src/PyRamanGUI.py'],
    pathex=[src_dir],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        'torch',
        'torchvision',
        'tensorflow',
        'pytest',
        'streamlit',
        'altair',
        'h5py',
        'beautifulsoup4',
        'selenium',
        'playwright',
        'behave',
        'black',
        'flake8',
        'isort',
        'mypy',
        'coverage',
        'pytest-cov',
        'pytest-html',
        'pytest-metadata',
        'pytest-ordering',
        'pytest-rerunfailures',
        'pytest-xdist',
        'seleniumbase',
        'fastapi',
        'uvicorn',
        'starlette',
        'weasyprint',
        'reportlab',
        'pyarrow',
        'plotly',
        'grpcio',
        'tensorboard',
        # cvxpy + native solvers: pulled in only by rampy >=0.6 (not the pinned
        # 0.4.9). They crash PyInstaller's isolated dependency-analysis subprocess
        # and are not used by PyRamanGUI -> exclude them.
        'cvxpy',
        'clarabel',
        'scs',
        'osqp',
        'highspy',
        'sparsediffpy',
        # no longer used after the OriginLab modules were removed
        'rampy',
        'sklearn',
        'scikit-learn',
        'pyqtgraph',
        'prettytable',
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='PyRaman_Pro',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,  # Set to False for GUI application
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,  # Add icon path here if you have one
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='PyRaman_Pro',
)
