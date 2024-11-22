import os
import subprocess
import ctypes
import numpy as np
import matplotlib.pyplot as plt


def test_construct_ishifted_transform_function():
    # Bind the C++ function
    transform_lib.construct_ishifted_transform_function.argtypes = [
        ctypes.POINTER(CmpxData),  # tf
        ctypes.c_int,             # n
        ctypes.c_float,           # wavelength
        ctypes.c_float,           # distance
        ctypes.c_float            # pixel_scale
    ]

    # Define parameters (typical for a space telescope simulation)
    n = 1024
    wavelength = 500e-9
    distance = 1000e-3
    pixel_scale = 10e-3 / (n / 2)

    # Allocate array for the output in Python
    tf_c = (CmpxData * (n * n))()

    # Call the C++ function
    transform_lib.construct_ishifted_transform_function(
        tf_c,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale)
    )

    # Convert C++ output to NumPy array
    tf_c_np = np.array([complex(val.real, val.imag) for val in tf_c]).reshape(n, n)

    # Python implementation for comparison
    delkx = 2.0 * np.pi / (pixel_scale * n)
    k = 2.0 * np.pi / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))
    scale = 1 / (n * n)
    tf_py = np.fft.ifftshift(scale * np.exp(1j * kz * distance))

    # Compute differences
    real_diff = np.abs(np.real(tf_py) - np.real(tf_c_np))
    imag_diff = np.abs(np.imag(tf_py) - np.imag(tf_c_np))

    relative_diff_imag = np.abs(imag_diff) / (np.abs(np.imag(tf_py)) + 1e-10)
    relative_diff_real = np.abs(real_diff) / (np.abs(np.real(tf_py)) + 1e-10)
    print(f"Max relative difference in imag component: {np.max(relative_diff_imag):.6e}")
    print(f"Max relative difference in real component: {np.max(relative_diff_real):.6e}")


    # Plot the results
    fig, axs = plt.subplots(3, 2, figsize=(15, 15))
    titles = [
        "Real Component (C++)", "Real Component (Python)",
        "Imag Component (C++)", "Imag Component (Python)",
        "Difference (Real)", "Difference (Imag)"
    ]
    data = [
        np.real(tf_c_np), np.real(tf_py),
        np.imag(tf_c_np), np.imag(tf_py),
        real_diff, imag_diff
    ]

    for ax, title, d in zip(axs.flat, titles, data):
        im = ax.imshow(d, extent=(0, n, 0, n), origin="lower", cmap="viridis")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.show()


def test_construct_transform_function():
    # Bind the C++ function
    transform_lib.construct_transform_function.argtypes = [
        ctypes.POINTER(CmpxData),  # tf
        ctypes.c_int,             # n
        ctypes.c_float,           # wavelength
        ctypes.c_float,           # distance
        ctypes.c_float            # pixel_scale
    ]

    # Define parameters (typical for a space telescope simulation)
    n = 1024
    wavelength = 500e-9
    distance = 1000e-3
    pixel_scale = 10e-3 / (n / 2)

    # Allocate array for the output in Python
    tf_c = (CmpxData * (n * n))()

    # Call the C++ function
    transform_lib.construct_transform_function(
        tf_c,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale)
    )

    # Convert C++ output to NumPy array
    tf_c_np = np.array([complex(val.real, val.imag) for val in tf_c]).reshape(n, n)

    # Python implementation for comparison
    delkx = 2.0 * np.pi / (pixel_scale * n)
    k = 2.0 * np.pi / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))
    scale = 1 / (n * n)
    tf_py = scale * np.exp(1j * kz * distance)

    # Compute differences
    real_diff = np.abs(np.real(tf_py) - np.real(tf_c_np))
    imag_diff = np.abs(np.imag(tf_py) - np.imag(tf_c_np))
    relative_diff_imag = np.abs(imag_diff) / (np.abs(np.imag(tf_py)) + 1e-10)
    relative_diff_real = np.abs(real_diff) / (np.abs(np.real(tf_py)) + 1e-10)
    print(f"Max relative difference in imag component: {np.max(relative_diff_imag):.6e}")
    print(f"Max relative difference in real component: {np.max(relative_diff_real):.6e}")

    # Plot the results
    fig, axs = plt.subplots(3, 2, figsize=(15, 15))
    titles = [
        "Real Component (C++)", "Real Component (Python)",
        "Imag Component (C++)", "Imag Component (Python)",
        "Difference (Real)", "Difference (Imag)"
    ]
    data = [
        np.real(tf_c_np), np.real(tf_py),
        np.imag(tf_c_np), np.imag(tf_py),
        real_diff, imag_diff
    ]

    for ax, title, d in zip(axs.flat, titles, data):
        im = ax.imshow(d, extent=(0, n, 0, n), origin="lower", cmap="viridis")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.show()





   # Bind the C++ function
    transform_lib.construct_transform_function.argtypes = [
        ctypes.POINTER(CmpxData),  # tf
        ctypes.c_int,             # n
        ctypes.c_float,           # wavelength
        ctypes.c_float,           # distance
        ctypes.c_float            # pixel_scale
    ]

    # Define parameters (typical for a space telescope simulation)
    n = 1024
    wavelength = 500e-9
    distance = 1000e-3
    pixel_scale = 10e-3 / (n / 2)

    # Allocate array for the output in Python
    tf_c = (CmpxData * (n * n))()

    # Call the C++ function
    transform_lib.construct_transform_function(
        tf_c,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale)
    )

    # Convert C++ output to NumPy array
    tf_c_np = np.array([complex(val.real, val.imag) for val in tf_c]).reshape(n, n)

    # Python implementation for comparison
    delkx = 2.0 * np.pi / (pixel_scale * n)
    k = 2.0 * np.pi / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))
    scale = 1 / (n * n)
    tf_py = scale * np.exp(1j * kz * distance)

    # Compute differences
    real_diff = np.abs(np.real(tf_py) - np.real(tf_c_np))
    imag_diff = np.abs(np.imag(tf_py) - np.imag(tf_c_np))
    relative_diff_imag = np.abs(imag_diff) / (np.abs(np.imag(tf_py)) + 1e-10)
    relative_diff_real = np.abs(real_diff) / (np.abs(np.real(tf_py)) + 1e-10)
    print(f"Max relative difference in imag component: {np.max(relative_diff_imag):.6e}")
    print(f"Max relative difference in real component: {np.max(relative_diff_real):.6e}")

    # Plot the results
    fig, axs = plt.subplots(3, 2, figsize=(15, 15))
    titles = [
        "Real Component (C++)", "Real Component (Python)",
        "Imag Component (C++)", "Imag Component (Python)",
        "Difference (Real)", "Difference (Imag)"
    ]
    data = [
        np.real(tf_c_np), np.real(tf_py),
        np.imag(tf_c_np), np.imag(tf_py),
        real_diff, imag_diff
    ]

    for ax, title, d in zip(axs.flat, titles, data):
        im = ax.imshow(d, extent=(0, n, 0, n), origin="lower", cmap="viridis")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.show()
    

def generate_star_gaussian(size, sigma, intensity, noise_stddev=0.1):
    # Create an empty image
    image = np.zeros((size, size))
    
    # Define the center of the image
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

    # Add random noise (Gaussian noise)
    noise = np.random.normal(0, noise_stddev, size=(size, size))
    image += noise

    # Clip the image to be in the valid range [0, 1]
    image = np.clip(image, 0, 1)

    return image


def test_angular_spectrum_propagation():

    # Bind the C++ function
    transform_lib.angular_spectrum_propagation.argtypes = [
        ctypes.POINTER(CmpxData),   # wavefront
        ctypes.c_int,               # n
        ctypes.c_float,             # wavelength
        ctypes.c_float,             # distance
        ctypes.c_float              # pixel_scale
    ]

    # Define parameters (typical for a space telescope simulation)
    n = 1024
    wavelength = 500e-9
    distance = 1000e-3
    pixel_scale = 10e-3 / (n / 2)

    # Generate Gaussian wavefront
    wavefront_c = (CmpxData * (n * n))()
    wavefront = generate_star_gaussian(n, 1, 1e5, 0.01) # generate_random_wavefront(n)
    for i in range(n * n):
        wavefront_c[i].real = wavefront.flat[i].real
        wavefront_c[i].imag = wavefront.flat[i].imag

    # Call the C++ function
    transform_lib.angular_spectrum_propagation(
        wavefront_c,
        n,
        ctypes.c_float(wavelength),
        ctypes.c_float(distance),
        ctypes.c_float(pixel_scale)
    )

    # Convert C++ output to NumPy array
    wavefront_c_np = np.array([complex(val.real, val.imag) for val in wavefront_c]).reshape(n, n)

    # Python implementation for comparison
    wf_as = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(wavefront)))  # Fourier transform of the wavefront
    delkx = 2.0 * np.pi / (pixel_scale * n)
    k = 2.0 * np.pi / wavelength
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))
    scale = 1
    tf_py = np.exp(1j * kz * distance)
    tf_py = scale * tf_py
    prop_wf = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(wf_as * tf_py)))  # Propagate the wavefront

    # Compute differences
    real_diff = np.abs(np.real(prop_wf) - np.real(wavefront_c_np))
    imag_diff = np.abs(np.imag(prop_wf) - np.imag(wavefront_c_np))

    # Calculate global min and max for real and imaginary components
    real_min = min(np.min(np.real(wavefront_c_np)), np.min(np.real(prop_wf)))
    real_max = max(np.max(np.real(wavefront_c_np)), np.max(np.real(prop_wf)))

    imag_min = min(np.min(np.imag(wavefront_c_np)), np.min(np.imag(prop_wf)))
    imag_max = max(np.max(np.imag(wavefront_c_np)), np.max(np.imag(prop_wf)))

    # Plot the results
    fig, axs = plt.subplots(3, 2, figsize=(15, 15))
    titles = [
        "Real Component (C++)",     "Imag Component (C++)",
        "Real Component (Python)",  "Imag Component (Python)",
        "Difference (Real)",        "Difference (Imag)"
    ]
    data = [
        np.real(wavefront_c_np), np.imag(wavefront_c_np), 
        np.real(prop_wf), np.imag(prop_wf),
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



if __name__ == "__main__":
    # Compile the C++ shared library
    cpp_files = "angular_spectrum_propagation.cpp fft.cpp"
    output_lib = "transform.so"
    compile_command = f"g++ -std=c++17 -shared -fPIC -o {output_lib} {cpp_files}"

    print("Compiling C++ library...")
    try:
        subprocess.run(compile_command, shell=True, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError as e:
        print(f"Compilation failed: {e}")
        exit(1)

    # Define the complex data type equivalent to cmpx_data_t
    class CmpxData(ctypes.Structure):
        _fields_ = [("real", ctypes.c_float), ("imag", ctypes.c_float)]

    # Load the shared library
    transform_lib = ctypes.CDLL(f'./{output_lib}')

    # test_construct_ishifted_transform_function()
    # test_construct_transform_function()
    test_angular_spectrum_propagation()
