"""
Generate synthetic Raman spectra for testing PyRamanGUI.

Creates spectra with:
- Realistic wavenumber range (100-3500 cm^-1)
- Lorentzian peaks at typical Raman positions
- Polynomial baseline (simulating fluorescence)
- Gaussian noise
- Optional cosmic spike
"""

import numpy as np
import os

def lorentzian(x, amplitude, center, width):
    """Lorentzian peak shape."""
    return amplitude * (width/2)**2 / ((x - center)**2 + (width/2)**2)

def gaussian(x, amplitude, center, width):
    """Gaussian peak shape."""
    return amplitude * np.exp(-4 * np.log(2) * ((x - center) / width)**2)

def generate_baseline(x, coeffs):
    """Generate polynomial baseline."""
    baseline = np.zeros_like(x)
    for i, c in enumerate(coeffs):
        baseline += c * (x / 1000)**i
    return baseline

def generate_spectrum(
    x,
    peaks,
    baseline_coeffs=(100, 50, 10),
    noise_level=5,
    add_cosmic_spike=False,
    spike_position=None
):
    """
    Generate a synthetic Raman spectrum.

    Parameters
    ----------
    x : array
        Wavenumber array
    peaks : list of tuples
        Each tuple: (amplitude, center, width, shape)
        shape: 'lorentzian' or 'gaussian'
    baseline_coeffs : tuple
        Polynomial coefficients for baseline
    noise_level : float
        Standard deviation of Gaussian noise
    add_cosmic_spike : bool
        Whether to add a cosmic spike
    spike_position : int or None
        Index for cosmic spike (random if None)

    Returns
    -------
    y : array
        Intensity values
    """
    # Start with baseline
    y = generate_baseline(x, baseline_coeffs)

    # Add peaks
    for amp, center, width, shape in peaks:
        if shape == 'lorentzian':
            y += lorentzian(x, amp, center, width)
        elif shape == 'gaussian':
            y += gaussian(x, amp, center, width)

    # Add noise
    y += np.random.normal(0, noise_level, len(x))

    # Add cosmic spike
    if add_cosmic_spike:
        if spike_position is None:
            spike_position = np.random.randint(100, len(x) - 100)
        y[spike_position] += 500 + np.random.normal(0, 50)

    return y

def save_spectrum(filename, x, y, header=None):
    """Save spectrum to text file."""
    if header:
        np.savetxt(filename, np.column_stack([x, y]),
                   header=header, delimiter='\t', fmt='%.6f')
    else:
        np.savetxt(filename, np.column_stack([x, y]),
                   delimiter='\t', fmt='%.6f')
    print(f"Saved: {filename}")

# Create output directory
os.makedirs(os.path.dirname(os.path.abspath(__file__)), exist_ok=True)
output_dir = os.path.dirname(os.path.abspath(__file__))

# Common wavenumber range
x = np.linspace(100, 3500, 3401)

# =============================================================================
# Spectrum 1: Silicon reference (single sharp peak at 520 cm^-1)
# =============================================================================
silicon_peaks = [
    (1000, 520, 5, 'lorentzian'),  # Main Si peak
]
y_silicon = generate_spectrum(x, silicon_peaks,
                               baseline_coeffs=(50, 5, 0),
                               noise_level=3)
save_spectrum(os.path.join(output_dir, 'silicon_reference.txt'), x, y_silicon,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 2: Graphene/Carbon (D, G, 2D peaks)
# =============================================================================
carbon_peaks = [
    (300, 1350, 40, 'lorentzian'),   # D band
    (800, 1580, 25, 'lorentzian'),   # G band
    (600, 2700, 50, 'lorentzian'),   # 2D band
]
y_carbon = generate_spectrum(x, carbon_peaks,
                              baseline_coeffs=(100, 30, 5),
                              noise_level=5)
save_spectrum(os.path.join(output_dir, 'carbon_graphene.txt'), x, y_carbon,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 3: Complex spectrum with multiple overlapping peaks
# =============================================================================
complex_peaks = [
    (200, 300, 20, 'lorentzian'),
    (400, 450, 30, 'lorentzian'),
    (350, 480, 25, 'lorentzian'),   # Overlapping with previous
    (500, 800, 40, 'gaussian'),
    (300, 1000, 35, 'lorentzian'),
    (250, 1050, 30, 'lorentzian'),  # Shoulder peak
    (600, 1500, 50, 'lorentzian'),
    (150, 2000, 60, 'gaussian'),
    (200, 2900, 80, 'gaussian'),    # C-H stretch region
    (180, 3000, 70, 'gaussian'),
]
y_complex = generate_spectrum(x, complex_peaks,
                               baseline_coeffs=(200, 80, 20, 2),
                               noise_level=8)
save_spectrum(os.path.join(output_dir, 'complex_spectrum.txt'), x, y_complex,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 4: High fluorescence background
# =============================================================================
fluorescence_peaks = [
    (150, 500, 30, 'lorentzian'),
    (200, 1000, 40, 'lorentzian'),
    (100, 1500, 35, 'lorentzian'),
]
y_fluorescence = generate_spectrum(x, fluorescence_peaks,
                                    baseline_coeffs=(500, 300, 100, 20),
                                    noise_level=10)
save_spectrum(os.path.join(output_dir, 'high_fluorescence.txt'), x, y_fluorescence,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 5: With cosmic spike
# =============================================================================
spike_peaks = [
    (400, 520, 10, 'lorentzian'),
    (300, 1000, 30, 'lorentzian'),
]
y_spike = generate_spectrum(x, spike_peaks,
                             baseline_coeffs=(80, 20, 5),
                             noise_level=5,
                             add_cosmic_spike=True,
                             spike_position=1500)
save_spectrum(os.path.join(output_dir, 'cosmic_spike.txt'), x, y_spike,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 6: Low SNR (noisy)
# =============================================================================
noisy_peaks = [
    (100, 500, 30, 'lorentzian'),
    (150, 1000, 40, 'lorentzian'),
    (80, 1500, 35, 'lorentzian'),
]
y_noisy = generate_spectrum(x, noisy_peaks,
                             baseline_coeffs=(100, 30, 10),
                             noise_level=30)
save_spectrum(os.path.join(output_dir, 'low_snr_noisy.txt'), x, y_noisy,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Spectrum 7: Clean spectrum (high SNR, for baseline testing)
# =============================================================================
clean_peaks = [
    (500, 400, 15, 'lorentzian'),
    (800, 800, 20, 'lorentzian'),
    (600, 1200, 25, 'lorentzian'),
    (400, 1600, 30, 'lorentzian'),
]
y_clean = generate_spectrum(x, clean_peaks,
                             baseline_coeffs=(150, 50, 15, 3),
                             noise_level=2)
save_spectrum(os.path.join(output_dir, 'clean_high_snr.txt'), x, y_clean,
              header='Wavenumber(cm-1)\tIntensity')

# =============================================================================
# Multiple spectra for batch processing / PCA testing
# =============================================================================
for i in range(5):
    # Vary peak intensities slightly
    batch_peaks = [
        (400 + np.random.normal(0, 30), 520, 10, 'lorentzian'),
        (300 + np.random.normal(0, 20), 1000, 30 + i*2, 'lorentzian'),
        (200 + i*30, 1500, 40, 'lorentzian'),  # This peak grows with i
    ]
    y_batch = generate_spectrum(x, batch_peaks,
                                 baseline_coeffs=(100 + i*10, 30, 5),
                                 noise_level=5)
    save_spectrum(os.path.join(output_dir, f'batch_spectrum_{i+1}.txt'), x, y_batch,
                  header='Wavenumber(cm-1)\tIntensity')

print(f"\nGenerated test spectra in: {output_dir}")
print("\nTest files created:")
print("  - silicon_reference.txt     : Simple single peak (520 cm-1)")
print("  - carbon_graphene.txt       : D, G, 2D bands")
print("  - complex_spectrum.txt      : Multiple overlapping peaks")
print("  - high_fluorescence.txt     : Strong baseline")
print("  - cosmic_spike.txt          : Contains cosmic ray spike")
print("  - low_snr_noisy.txt         : High noise level")
print("  - clean_high_snr.txt        : Clean spectrum for baseline testing")
print("  - batch_spectrum_1-5.txt    : Series for PCA/batch testing")
