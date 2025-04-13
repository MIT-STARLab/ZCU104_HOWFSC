
import numpy as np

def fft(x):
  """
  A recursive implementation of 
  the 1D Cooley-Tukey FFT, the 
  input should have a length of 
  power of 2. 
  """
  N = len(x)
  
  if N == 1:
      return x
  
  X_even = fft(x[::2])
  X_odd = fft(x[1::2])

  factor = np.exp(-2j*np.pi*np.arange(N)/ N)
  
  X = np.concatenate([X_even+factor[:int(N/2)]*X_odd, X_even+factor[int(N/2):]*X_odd])
  
  return X


def fft_numpy(x):
   return np.fft.fft(x)


def fft2d(x):
    """
    A 2D FFT is implemented by doing a 1D FFT on all the rows in the input array 
    then doing a 1DFFT on the cols in the resulting matrix
    """
    m,n = len(x), len(x[0])
    output = np.zeros((m, n), dtype=complex)

    # first do a 1D FFT on all the cols
    for row_index in range(m):
        output[row_index] = fft(x[row_index])

    for col_index in range(n):
        output[:,col_index] = fft(output[:,col_index])

    return output


def fft2d_numpy(x):
    return np.fft.fft2(x)


def compare_arrays(arr1, arr2, tolerance):

    # quick fails assertions
    assert len(arr1) == len(arr2), 'arrays must have the same length'
    if len(arr1) == 0:
        return True

    m = len(arr1[0])
    for i in range(len(arr1)):
        assert len(arr1[i]) == m
        assert len(arr2[i]) == m

    m, n = len(arr1), len(arr1[0])

    # actual checking:
    passed = True
    for i in range(m):
        for j in range(n):
            if np.abs(arr1[i][j].real-arr2[i][j].real) > tolerance or np.abs(arr1[i][j].imag-arr2[i][j].imag) > tolerance:
                print(f"Error at location {(i,j)}, first value is {arr1[i][j]}, second value is {arr2[i][j]}")
                passed = False
        
    return passed


def testing_fft2d():

    trials = 1
    m = 2**10
    n = m*2
    tolerance = 1e-10

    for i in range(trials):
        x = np.random.rand(m, n) + 1j * 0
        numpy_output = np.fft.fft2(x)
        my_outpuy = fft2d(x)

        print("Trial 1 results:")
        print(compare_arrays(numpy_output, my_outpuy, tolerance))


testing_fft2d()