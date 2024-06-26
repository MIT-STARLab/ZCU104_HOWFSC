import numpy as np
import matplotlib.pyplot as plt

# Parameters
FFT_LENGTH = 2**10
sampling_interval = 0.01

# Time array
t = np.arange(0, sampling_interval * FFT_LENGTH, sampling_interval)

# Generate signal
freq = 1.
x_real = 3 * np.sin(2 * np.pi * freq * t)

freq = 4
x_real += np.sin(2 * np.pi * freq * t)

freq = 7   
x_real += 0.5 * np.sin(2 * np.pi * freq * t)

# Set imaginary part of x to 0
x_imag = np.zeros_like(x_real)

# Combine x_real and x_imag into x
x = x_real + 1j * x_imag

# Perform FFT
y = np.fft.fft(x)

# Save x and y to CSV
data = np.column_stack((x_real, x_imag, y.real, y.imag))  # Combine x_real, x_imag, real part of y, and imaginary part of y
np.savetxt('ftt_tb_data_amd.csv', data, delimiter=',', header='x_real, x_imag, y_real, y_imag', comments='')

# Plotting (optional)
# plt.figure(figsize=(10, 6))
# plt.plot(t, x_real, label='Signal x_real')
# plt.legend()
# plt.xlabel('Time')
# plt.ylabel('Amplitude')
# plt.title('Signal x_real vs Time')
# plt.grid(True)
# plt.show()