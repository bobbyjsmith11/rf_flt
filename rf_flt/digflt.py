#!/usr/bin/env python
"""
===========================
digflt.py
===========================
:Description:
    various calculations and plotting involving digital filters

"""


import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def load_coeffs_from_file(fpath):
    """
    """
    fo = open(fpath, 'r')
    s = fo.read()
    coeffs = []
    for c in s.split(","):
        coeffs.append(int(c))
    return coeffs

def plot_fir(coeffs=[1, -1], sr=44100, normalized=True, numbits=16):
    """
    """
    plt.ion()
    coeffs = np.array(coeffs)
    if normalized:
        coeffs = coeffs/(2**(numbits-1))
    w, h = signal.freqz(b=coeffs, a=1)
    x = w * sr * 1.0 / (2 * np.pi)
    y = 20 * np.log10(abs(h))
    fig = plt.figure(figsize=(10,5))
    # plt.semilogx(x, y)
    plt.plot(x, y)
    plt.ylabel('Amplitude [dB]')
    plt.xlabel('Frequency [Hz]')
    title = 'Frequency response\n{:0.6f} MSPS'.format(sr/1e6)
    plt.title(title)
    plt.grid(which='both', linestyle='-', color='grey')
    plt.xlim(0, sr/2)
    plt.show()
    return fig
