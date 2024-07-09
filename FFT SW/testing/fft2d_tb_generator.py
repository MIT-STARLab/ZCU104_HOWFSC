
import numpy as np

m = 2**10
n = m
input_array = np.random.rand(m, n) + 1j * 0
numpy_output = np.fft.fft2(input_array)

input_array = input_array.flatten('C')
numpy_output = numpy_output.flatten('C')

data = np.column_stack((input_array.real, input_array.imag, numpy_output.real, numpy_output.imag))
np.savetxt('fft2d_tb_data.csv', data, delimiter=',', header='x_real, x_imag, y_real, y_imag', comments='')
