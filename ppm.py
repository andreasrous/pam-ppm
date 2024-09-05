import sys
import matplotlib.pyplot as plt
from functions import *

def plot_name_to_ppm_waveform(name, M):
    # Βήμα 1: μετατροπή του ονόματος σε σειρά από 8 bits βάσει του κώδικα ASCII
    bits = name_to_binary(name)

    # Βήμα 2: μετατροπή των bits σε σύμβολα M-PPM βάσει του κώδικα Gray μήκους log2M
    fmap = pam_gray_forward_map(M)
    symbols = bits_to_symbols(bits, fmap)

    # Βήμα 3: δημιουργία κυματομορφής PPM
    TS = 1e-9  # Διάρκεια συμβόλου ώστε Rb = 1Gb/s (ισχύει Rb = 1/TS)
    t, x = ppm_waveform(symbols, TS, M, fmap)
    plot_ppm(t, x, M, bits, TS, dy = 1, title = f"M = {M}")

if len(sys.argv) < 2:
    print("Usage: ppm.py <full_name>")
else:
    full_name = ' '.join(sys.argv[1:])
    for M in [2, 4, 8, 16]:
        plot_name_to_ppm_waveform(full_name, M)

# Show all figures
plt.show()
