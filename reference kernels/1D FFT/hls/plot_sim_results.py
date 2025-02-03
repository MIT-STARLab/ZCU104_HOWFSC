import matplotlib.pyplot as plt
import numpy as np

def read_data(filename):
    real = []
    imag = []
    with open(filename, 'r') as file:
        for line in file:
            r, i = map(float, line.strip().split())
            real.append(r)
            imag.append(i)
    return np.array(real), np.array(imag)

# Plot input, real, and imaginary parts in separate subplots
plt.figure(figsize=(12, 15))

# Read the data
input_real, input_imag = read_data('fft_krnl/hls/sim/wrapc/software_input.txt')
software_real, software_imag = read_data('fft_krnl/hls/sim/wrapc/software_output.txt')
hardware_real, hardware_imag = read_data('fft_krnl/hls/sim/wrapc/hardware_output.txt')

# Input Real Part
plt.subplot(5, 1, 1)
plt.plot(input_real, label='Input Real', color='black', alpha=0.5)
plt.legend()
plt.ylabel('Amplitude')

# Input Imaginary Part#
# plt.subplot(5, 1, 2)
# plt.plot(input_imag, label='Input Imaginary', color='black', linestyle='--', alpha=0.5)
# plt.legend()
# plt.ylabel('Amplitude')

# Real part of Software FFT
plt.subplot(5, 1, 2)
plt.plot(software_real/2048, label='Software FFT Real', color='blue', alpha=0.5)
plt.legend()
plt.ylabel('Amplitude')

# Imaginary part of Software FFT
plt.subplot(5, 1, 3)
plt.plot(software_imag/2048, label='Software FFT Imaginary', color='blue', alpha=0.5)
plt.legend()
plt.ylabel('Amplitude')

# Real part of Hardware FFT
plt.subplot(5, 1, 4)
plt.plot(hardware_real/2048, label='Hardware FFT Real', color='red', alpha=0.5)
plt.legend()
plt.ylabel('Amplitude')

# Imaginary part of Hardware FFT
plt.subplot(5, 1, 5)
plt.plot(hardware_imag/2048, label='Hardware FFT Imaginary', color='blue', alpha=0.5)
plt.legend()
plt.ylabel('Amplitude')
plt.xlabel('Sample Index')

plt.tight_layout()
plt.show()