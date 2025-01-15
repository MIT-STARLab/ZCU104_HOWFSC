"""
 " MIT STAR Lab
 " M.Subhi Abo Rdan (msubhi_a@mit.edu)
 " Last modified on January 15, 2024
 " Pure software implementations of Angular Spectrum Propagation Method + 
 " Scripts to test C++ implementations
"""

import os
import subprocess
import ctypes
import numpy as np
import matplotlib.pyplot as plt


########### Software Implmentations ###########

def construct_propagation_matrix(n, pixel_scale, wavelength, distance, scale = False, shift = False):
    """
    TODO: docs
    """
    delkx = 2.0 * np.pi / (pixel_scale * n)
    k = 2.0 * np.pi / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))
    transfer_function = np.exp(1j * kz * distance)

    if scale:
        transfer_function = transfer_function * (1/(n**2))

    if shift:
        transfer_function = np.fft.ifftshift(transfer_function)

    return transfer_function



def angular_spectrum_progagation(wavefront, n, pixel_scale, wavelength, distance):
    """
    TODO: docs
    """

    transfer_function = construct_propagation_matrix(n, pixel_scale, wavelength, distance, scale = True, shift = True)

    wf_as = np.fft.fft2(np.fft.fftshift(wavefront), norm = "backward")
    wf_as = wf_as * transfer_function
    wf_as = np.fft.fftshift(np.fft.ifft2(wf_as, norm = "forward"))  # Scaling in Transfer Function

    return wf_as



def angular_spectrum_progagation_float(wavefront, n, pixel_scale, wavelength, distance):
    """
    TODO: docs
    """
    # Ensure inputs are float32
    wavefront = wavefront.astype(np.float32)
    pixel_scale = np.float32(pixel_scale)
    wavelength = np.float32(wavelength)
    distance = np.float32(distance)

    # Precompute constants in float32
    delkx = np.float32(2.0 * np.pi) / (pixel_scale * n)
    k = np.float32(2.0 * np.pi) / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)], dtype=np.float32)

    # Create float32 meshgrids
    kx, ky = np.meshgrid(kxy, kxy, indexing="ij")
    kx = kx.astype(np.float32)
    ky = ky.astype(np.float32)

    # Compute kz with float32 precision
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j)).astype(np.complex64)
    scale = np.float32(1.0) / (n**2)
    tf_py = scale * np.exp(1j * kz * distance).astype(np.complex64)

    # Perform Fourier operations with single precision
    wf_as = np.fft.fft2(np.fft.fftshift(wavefront)).astype(np.complex64)
    wf_as = np.fft.ifftshift(wf_as)  # Fourier transform of the wavefront
    wf_as = wf_as * tf_py
    wf_as = np.fft.ifftshift(wf_as)

    # Use inverse Fourier transform
    wf_as = np.fft.fftshift(np.fft.ifft2(wf_as, norm="forward")).astype(np.complex64)

    return wf_as



def generate_star_gaussian(size, sigma, intensity, noise_stddev=0.1, noise = True):
    """
    TODO: docs
    """
    image = np.zeros((size, size))
    center_x, center_y = size // 2, size // 2

    # Generate Gaussian pattern for the first star
    x_grid, y_grid = np.meshgrid(np.arange(size), np.arange(size))
    distance = np.sqrt((x_grid - center_x)**2 + (y_grid - center_y)**2)
    gaussian = intensity * np.exp(-distance**2 / (2 * sigma**2))
    image += gaussian

    # Generate Gaussian pattern for the second star with smaller intensity
    x_grid, y_grid = np.meshgrid(np.arange(size), np.arange(size))
    distance = np.sqrt((x_grid - center_x + size // 8)**2 + (y_grid - center_y + size // 8)**2)
    gaussian = (intensity / 8) * np.exp(-distance**2 / (2 * (sigma/2)**2))
    image += gaussian

    if noise:
        # Add random noise (Gaussian noise)
        noise_data = np.random.normal(0, noise_stddev, size=(size, size))
        image += noise_data

    return image



########### Testing Utils ###########

def compute_error(arr1, arr2):
    """
    computes and prints the average relative error in both real and image componenets
    """
    real_arr1 = np.real(arr1)
    imag_arr1 = np.imag(arr1)

    real_arr2 = np.real(arr2)
    imag_arr2 = np.imag(arr2)

    # Compute absolute and relative errors for real and imaginary parts
    relative_error_real = np.abs(real_arr1 - real_arr2) / (np.abs(real_arr1) + 1e-12)  # Add epsilon to avoid division by zero
    relative_error_imag = np.abs(imag_arr1 - imag_arr2) / (np.abs(imag_arr1) + 1e-12)

    # Compute average relative errors
    average_relative_error_real = np.mean(relative_error_real)
    average_relative_error_imag = np.mean(relative_error_imag)

    return average_relative_error_real, average_relative_error_imag



def comparisonPlot(arr1, arr2):
    """
    TODO: docs
    """
    
    # Compute differences
    real_diff = np.abs(np.real(arr1) - np.real(arr2))
    imag_diff = np.abs(np.imag(arr1) - np.imag(arr2))
    
    # Calculate global min and max for real and imaginary components
    real_min = min(np.min(np.real(arr2)), np.min(np.real(arr1)))
    real_max = max(np.max(np.real(arr2)), np.max(np.real(arr1)))

    imag_min = min(np.min(np.imag(arr2)), np.min(np.imag(arr1)))
    imag_max = max(np.max(np.imag(arr2)), np.max(np.imag(arr1)))

    # Plot the results
    fig, axs = plt.subplots(3, 2, figsize=(15, 15))
    titles = [
        "Real Component (C++)",     "Imag Component (C++)",
        "Real Component (Python)",  "Imag Component (Python)",
        "Difference (Real)",        "Difference (Imag)"
    ]
    data = [
        np.real(arr2), np.imag(arr2), 
        np.real(arr1), np.imag(arr1),
        real_diff, imag_diff
    ]

    # Create consistent colormap scaling
    for i, (ax, title, d) in enumerate(zip(axs.flat, titles, data)):
        if "Difference" in title:
            cmap = "coolwarm"  # Different colormap for differences
            vmin, vmax = np.min(d), np.max(d)  # Use min/max for the difference
        else:
            cmap = "viridis"  # Consistent colormap for real and imaginary components
            if "Real" in title:
                vmin, vmax = real_min, real_max  # Same scale for real components
            else:
                vmin, vmax = imag_min, imag_max  # Same scale for imaginary components

        im = ax.imshow(d, extent=(0, n, 0, n), origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title)
        fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.show()



########### Tests ###########

def test_angular_spectrum_propagation_float(n, wavelength, distance, pixel_scale):

    # Bind the C++ function
    transform_lib.angular_spectrum_propagation_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),   # wavefront
        ctypes.c_int,               # n
        ctypes.c_float,             # wavelength
        ctypes.c_float,             # distance
        ctypes.c_float              # pixel_scale
    ]

    # Generate Gaussian wavefront
    wavefront_py = generate_star_gaussian(n, 100 , 1e5, 0.01, False)

    wavefront_cpp = (CmpxData * (n * n))()
    for i in range(n * n):
        wavefront_cpp[i].real = wavefront_py.flat[i].real
        wavefront_cpp[i].imag = wavefront_py.flat[i].imag

    # Call the C++ function and convert output to NumPy array
    transform_lib.angular_spectrum_propagation_float_wrapper(
        wavefront_cpp,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale)
    )
    wavefront_cpp = np.array([complex(val.real, val.imag) for val in wavefront_cpp]).reshape(n, n)
    wavefront_py = angular_spectrum_progagation(wavefront_py, n, pixel_scale, wavelength, distance)

    comparisonPlot(wavefront_py, wavefront_cpp)



def test_construct_ishifted_transform_function_float(n, wavelength, distance, pixel_scale):

    # Bind the C++ function
    transform_lib.construct_ishifted_transform_function_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),  # tf
        ctypes.c_int,              # n
        ctypes.c_float,            # wavelength
        ctypes.c_float,            # distance
        ctypes.c_float,            # pixel_scale
        ctypes.c_bool              # mod
    ]

    tf_py = construct_propagation_matrix(n, pixel_scale, wavelength, distance, scale =  True, shift =  True)

    tf_cpp = (CmpxData * (n * n))()
    transform_lib.construct_ishifted_transform_function_float_wrapper(
        tf_cpp,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale),
        ctypes.c_bool(True)
    )
    tf_cpp = np.array([complex(val.real, val.imag) for val in tf_cpp]).reshape(n, n)

    comparisonPlot(tf_py, tf_cpp)



def test_construct_transform_function_float(n, wavelength, distance, pixel_scale):

    # Bind the C++ function
    transform_lib.construct_transform_function_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),  # tf
        ctypes.c_int,              # n
        ctypes.c_float,            # wavelength
        ctypes.c_float,            # distance
        ctypes.c_float,            # pixel_scale
        ctypes.c_bool              # mod
    ]

    tf_py = construct_propagation_matrix(n, pixel_scale, wavelength, distance, scale =  True, shift =  False)

    tf_cpp = (CmpxData * (n * n))()
    transform_lib.construct_transform_function_float_wrapper(
        tf_cpp,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale),
        ctypes.c_bool(True)
    )
    tf_cpp = np.array([complex(val.real, val.imag) for val in tf_cpp]).reshape(n, n)

    comparisonPlot(tf_py, tf_cpp)



def test_generate_star_gaussian_floats():

    transform_lib.generate_star_gaussian_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),   # wavefront
        ctypes.c_int,               # size
        ctypes.c_float,             # sigma
        ctypes.c_float,             # intensity
        ctypes.c_float,             # noise_stddev
        ctypes.c_bool,              # noise
    ]

    # Define parameters (typical for a space telescope simulation)
    n = 1024
    sigma = 100
    intensity = 1e5
    noise_stddev = 0.01
    noise = False

    wavefront_c = (CmpxData * (n * n))()
    transform_lib.generate_star_gaussian_float_wrapper(
        wavefront_c,
        n,
        ctypes.c_float(sigma),
        ctypes.c_float(intensity),
        ctypes.c_float(noise_stddev),
        ctypes.c_bool(noise)
    )
    wavefront_c_real = np.array([wavefront_c[i].real for i in range(n * n)]).reshape(n, n)
    wavefront_python = generate_star_gaussian(n, sigma, intensity, noise_stddev, noise)


    # Compute the difference
    wavefront_diff = wavefront_c_real - wavefront_python

    # Plot the results
    plt.figure(figsize=(15, 5))

    plt.subplot(1, 3, 1)
    plt.title("C Function Output (Real)")
    plt.imshow(wavefront_c_real, cmap="gray")
    plt.colorbar()

    plt.subplot(1, 3, 2)
    plt.title("Python Function Output")
    plt.imshow(wavefront_python, cmap="gray")
    plt.colorbar()

    plt.subplot(1, 3, 3)
    plt.title("Difference (Real Component)")
    plt.imshow(wavefront_diff, cmap="gray")
    plt.colorbar()

    plt.tight_layout()
    plt.show()



def test_double_floats_difference_python(n, wavelength, distance, pixel_scale):
    wavefront = generate_star_gaussian(n, 100 , 1e5, 0.01, False)

    # Run using two different functions
    double_results = angular_spectrum_progagation(wavefront, n, pixel_scale, wavelength, distance)
    float_results  = angular_spectrum_progagation_float(wavefront, n, pixel_scale, wavelength, distance)

    # compute average relative error
    average_relative_error_real, average_relative_error_imag = compute_error(double_results, float_results)
    print(f"Average Relative Error (Real Component): {average_relative_error_real:.6f}")
    print(f"Average Relative Error (Imag Component): {average_relative_error_imag:.6f}")

    comparisonPlot(double_results, float_results)




if __name__ == "__main__":

    ## compile the C++ source code and load the shared library
    cpp_files = "angular_spectrum_propagation.cpp fft.cpp"
    output_lib = "angular_spectrum_propagation.so"
    compile_command = f"g++ -DPYTHON_TESTING_FLOATS -std=c++17 -shared -fPIC -o {output_lib} {cpp_files}"

    print("Compiling C++ library...")
    try:
        subprocess.run(compile_command, shell=True, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed: {e}")
        exit(1)

    transform_lib = ctypes.CDLL(f'./{output_lib}')


    # Define the complex data type equivalent to complex<float>
    class CmpxData(ctypes.Structure):
        _fields_ = [("real", ctypes.c_float), ("imag", ctypes.c_float)]

    # Define testing parameters (typical for a space telescope simulation)
    n = 1024
    wavelength = 500e-9
    distance = 1000e-3
    pixel_scale = 10e-3 / (n / 2)

    ## Tests
    test_generate_star_gaussian_floats()
    test_construct_ishifted_transform_function_float(n, wavelength, distance, pixel_scale)
    test_construct_transform_function_float(n, wavelength, distance, pixel_scale)
    test_angular_spectrum_propagation_float(n, wavelength, distance, pixel_scale)
    test_double_floats_difference_python(n, wavelength, distance, pixel_scale)
