

import numpy as np
import matplotlib.pyplot as plt


n = 512
m = 128
file1 = "./file1.csv"
file2 = "./file2.csv"

def generate_fft_1d_testing_file(n, file1):
    # Parameters
    sampling_interval = 0.01

    # Time array
    t = np.arange(0, sampling_interval * n, sampling_interval)

    # Generate signal
    freq = 0.97
    x_real = 203 * np.sin(2 * np.pi * freq * t)

    freq = 10.44
    x_real += 124 * np.sin(2 * np.pi * freq * t)

    freq = 50.44
    x_real += 34 * np.sin(2 * np.pi * freq * t)

    # Set imaginary part of x to 0
    x_imag = np.zeros_like(x_real)

    # Combine x_real and x_imag into x
    x = x_real + 1j * x_imag

    # Perform FFT
    y = np.fft.fft(x)

    # Save x and y to CSV
    data = np.column_stack((x_real, x_imag, y.real, y.imag))  # Combine x_real, x_imag, real part of y, and imaginary part of y
    np.savetxt(file1, data, delimiter=',', header='x_real, x_imag, y_real, y_imag', comments='')



def generate_fft_2d_testing_file(m, n, file2):
    input_array = np.random.rand(m, n) + 1j * 0
    numpy_output = np.fft.fft2(input_array)

    input_array = input_array.flatten('C')
    numpy_output = numpy_output.flatten('C')

    data = np.column_stack((input_array.real, input_array.imag, numpy_output.real, numpy_output.imag))
    np.savetxt(file2, data, delimiter=',', header='x_real, x_imag, y_real, y_imag', comments='')


generate_fft_1d_testing_file(n, file1)
generate_fft_2d_testing_file(m, n, file2)
