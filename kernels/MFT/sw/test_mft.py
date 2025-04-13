"""
 " MIT STAR Lab
 " M.Subhi Abo Rdan (msubhi_a@mit.edu)
 " Last modified on Feb 23, 2025
 " Scripts to test MFT C++ implementations
"""

import subprocess
import ctypes
import numpy as np
import matplotlib.pyplot as plt



def mft_forward(wavefront, npix, npsf, psf_pixelscale_lamD, 
                convention='-', pp_centering='odd', fp_centering='odd'):
    """
    Code by Kian Milani
    """

    N = wavefront.shape[0]
    dx = 1.0 / npix
    if pp_centering=='even':
        Xs = (np.arange(N, dtype=float) - (N / 2) + 1/2) * dx
    elif pp_centering=='odd':
        Xs = (np.arange(N, dtype=float) - (N / 2)) * dx

    du = psf_pixelscale_lamD
    if fp_centering=='odd':
        Us = (np.arange(npsf, dtype=float) - npsf / 2) * du
    elif fp_centering=='even':
        Us = (np.arange(npsf, dtype=float) - npsf / 2 + 1/2) * du

    xu = np.outer(Us, Xs)
    vy = np.outer(Xs, Us)

    if convention=='-':
        My = np.exp(-1j*2*np.pi*vy) 
        Mx = np.exp(-1j*2*np.pi*xu)
    else:
        My = np.exp(1j*2*np.pi*vy) 
        Mx = np.exp(1j*2*np.pi*xu)

    norm_coeff = psf_pixelscale_lamD/npix


    return Mx@wavefront@My * norm_coeff


def mft_reverse(fpwf, psf_pixelscale_lamD, npix, N, convention='+', pp_centering='odd', fp_centering='odd'):
    """
    Code by Kian Milani
    """

    npsf = fpwf.shape[0]
    du = psf_pixelscale_lamD
    if fp_centering=='odd':
        Us = (np.arange(npsf, dtype=float) - npsf / 2) * du
    elif fp_centering=='even':
        Us = (np.arange(npsf, dtype=float) - npsf / 2 + 1/2) * du

    dx = 1.0 / npix
    if pp_centering=='even':
        Xs = (np.arange(N, dtype=float) - (N / 2) + 1/2) * dx
    elif pp_centering=='odd':
        Xs = (np.arange(N, dtype=float) - (N / 2)) * dx

    ux = np.outer(Xs, Us)
    yv = np.outer(Us, Xs)

    if convention=='+':
        My = np.exp(1j*2*np.pi*yv) 
        Mx = np.exp(1j*2*np.pi*ux)
    else:
        My = np.exp(-1j*2*np.pi*yv) 
        Mx = np.exp(-1j*2*np.pi*ux)

    norm_coeff = psf_pixelscale_lamD/npix 

    return Mx@fpwf@My * norm_coeff


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
    n = arr1.shape[0]
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




###################################          Tests              ##########################################

def test_mat_mul_float(m = 512, n = 128, p = 64):
    print(f"Testing complex matrix multiplication with random data using m = {m}, n={n}, p={p}")

    # Bind the C++ function
    # void complex_mat_mul_float_wrapper(std::complex<float> *C, std::complex<float> *A, std::complex<float> *B, int M, int N, int P);
    transform_lib.complex_mat_mul_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),   # C
        ctypes.POINTER(CmpxData),   # A
        ctypes.POINTER(CmpxData),   # B

        ctypes.c_int,              # M
        ctypes.c_int,              # N
        ctypes.c_int               # P
    ]


    A_py = np.random.rand(m, n).astype(float) + 1j * np.random.rand(m, n).astype(float)
    B_py = np.random.rand(n, p).astype(float) + 1j * np.random.rand(n, p).astype(float)
    C_py = np.dot(A_py, B_py)

    
    A_cpp = (CmpxData * (m * n))()
    B_cpp = (CmpxData * (n * p))()
    C_cpp = (CmpxData * (m * p))()

    # Copy Python complex arrays into C++ structures
    for i in range(m):
        for j in range(n):
            index = i*n + j
            A_cpp[index].real = A_py[i][j].real
            A_cpp[index].imag = A_py[i][j].imag

    for i in range(n):
        for j in range(p):
            index = i*p + j
            B_cpp[index].real = B_py[i][j].real
            B_cpp[index].imag = B_py[i][j].imag


    transform_lib.complex_mat_mul_float_wrapper(
        C_cpp,
        A_cpp,
        B_cpp,
        ctypes.c_int(m),
        ctypes.c_int(n),
        ctypes.c_int(p)
    )

    C_cpp = np.array([complex(val.real, val.imag) for val in C_cpp]).reshape(m, p)

    average_relative_error_real, average_relative_error_imag = compute_error(C_cpp, C_py)
    print(f"average_relative_error_real = {average_relative_error_real}")
    print(f"average_relative_error_imag = {average_relative_error_imag}")




def test_mft_forward_float(psf_pixelscale_lamD=0.1, npix=512, npsf=128, plot=False):
    print(f"Testing MFT Forward on gaussian input data with  psf_pixelscale_lamD = {psf_pixelscale_lamD}, npix = {npix}, npsf={npsf}")
    # Bind the C++ function
    # void mft_forward_float_wrapper(std::complex<float> *wavefront, std::complex<float> *wavefront_output, float psf_pixelscale_lamD, int npsf, int npix)
    transform_lib.mft_forward_float_wrapper.argtypes = [
        ctypes.POINTER(CmpxData),  # wavefront
        ctypes.POINTER(CmpxData),  # wavefront_output
        ctypes.c_float,            # psf_pixelscale_lamD
        ctypes.c_int,              # npsf
        ctypes.c_int               # npix
    ]

    # Generate the input wavefront
    wavefront_py = generate_star_gaussian(npix, 100, 1e5, 0.01, False)
    output_py = mft_forward(wavefront_py, npix, npsf, psf_pixelscale_lamD)

    # Allocate memory for C++ input
    wavefront_cpp = (CmpxData * (npix * npix))()
    for i in range(npix):
        for j in range(npix):
            idx = i * npix + j
            wavefront_cpp[idx].real = wavefront_py[i, j].real
            wavefront_cpp[idx].imag = wavefront_py[i, j].imag

    # Allocate memory for C++ output
    wavefront_output_cpp = (CmpxData * (npsf * npsf))()

    # Call the C++ function
    transform_lib.mft_forward_float_wrapper(
        wavefront_cpp,
        wavefront_output_cpp,
        ctypes.c_float(psf_pixelscale_lamD),
        ctypes.c_int(npsf),
        ctypes.c_int(npix)
    )

    # Convert C++ output to NumPy array
    output_cpp_np = np.array([complex(val.real, val.imag) for val in wavefront_output_cpp]).reshape(npsf, npsf)

    # Compute error and visualize
    average_relative_error_real, average_relative_error_imag = compute_error(output_cpp_np, output_py)
    print(f"average_relative_error_real = {average_relative_error_real}")
    print(f"average_relative_error_imag = {average_relative_error_imag}")

    if plot:
        comparisonPlot(output_py, output_cpp_np)






if __name__ == "__main__":

    ## compile the C++ source code and load the shared library
    cpp_files = "mft.cpp"
    output_lib = "mft.so"
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
    psf_pixelscale_lamD = 0.1
    npix = 2000
    npsf = 1000

    test_mat_mul_float()
    test_mft_forward_float(plot=True)
    print("Tests Done")
