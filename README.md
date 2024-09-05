# PAM and PPM Simulation in Python

This repository contains two Python scripts that simulate Pulse Amplitude Modulation (PAM) and Pulse Position Modulation (PPM) systems. The project was part of a university assignment for telecommunications systems.

## Files

1. **pam.py:** Contains functions to calculate the Signal-to-Noise Ratio (SNR) for PAM systems and validate the results.
2. **ppm.py:** Simulates the waveform generation for Pulse Position Modulation (PPM) and plots the signal based on a user's name.
3. **functions.py:** Includes utility functions shared between both `pam.py` and `ppm.py`, such as conversion functions and waveform plotting.

## Overview

### PAM (Pulse Amplitude Modulation)

- **PAM SNR Calculation:** The function `calculate_SNRbdB` computes the SNR for a given modulation order `M` and bit error probability `Pb`. This is validated using the `validate` function.
- **Plotting SNR:** The `plot_SNRbdB` function visualizes how SNR changes with different modulation orders `M` for a given `Pb`.

### PPM (Pulse Position Modulation)

- **PPM Waveform Generation:** The script `ppm.py` takes a name as input and generates PPM waveforms using the Gray code mapping technique. The waveform is plotted for different values of `M` (2, 4, 8, 16).

## Usage

### PAM Simulation

1. **Run the PAM SNR validation:**

   ```
   python pam.py
   ```

   This will calculate and validate the SNR for different modulation orders and bit error probabilities.

2. **Plot SNR vs Modulation Order:** Modify the `plot_SNRbdB` function's input to experiment with different values of `Pb`.

### PPM Simulation

1. **Generate PPM Waveforms:**

   ```
   python ppm.py <full_name>
   ```

   Replace `<full_name>` with the name you want to encode and generate the waveform for. The script will create and display the PPM waveform for modulation orders `M = 2, 4, 8, 16`.

## Dependencies

- **NumPy:** For mathematical operations.
- **Matplotlib:** For plotting the waveforms and graphs.
- **SciPy:** For the Q-function and error calculations.

You can install the required libraries with:

```
pip install numpy matplotlib scipy
```

## Example

To run the PPM simulation for the name "John Doe" with different modulation orders:

```
python ppm.py John Doe
```
