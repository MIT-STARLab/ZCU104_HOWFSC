import numpy as np

"""
Code by Kian Milani

Computes the forward Matrix Fourier Transform (MFT) of a given wavefront.

The MFT provides better control over output sampling compared to the standard Fast Fourier Transform (FFT),
making it useful in optical modeling and wavefront propagation simulations.

Parameters:
-----------
wavefront : np.ndarray (complex, shape: NxN)
    The input complex wavefront, typically representing an optical field.

npix : int
    The number of pixels in the input wavefront array.

npsf : int
    The number of pixels in the output focal plane (PSF sampling).

psf_pixelscale_lamD : float
    The pixel scale in the focal plane, typically in units of λ/D (lambda over diameter).

convention : str, optional (default='-')
    Defines the Fourier transform sign convention:
    - `'-'`: Uses `exp(-j 2π x)`, which is the standard convention for forward Fourier transforms.
    - `'+'`: Uses `exp(+j 2π x)`, which may be useful in specific cases.

pp_centering : str, optional (default='odd')
    Determines the centering convention for the pupil plane (`Xs`):
    - `'even'`: Centers with a half-pixel offset.
    - `'odd'`: Centers symmetrically.

fp_centering : str, optional (default='odd')
    Determines the centering convention for the focal plane (`Us`):
    - `'even'`: Centers with a half-pixel offset.
    - `'odd'`: Centers symmetrically.
"""


def mft_forward(wavefront, npix, npsf, psf_pixelscale_lamD, 
                convention='-', pp_centering='odd', fp_centering='odd'):

    N = wavefront.shape[0]
    dx = 1.0 / npix
    if pp_centering=='even':
        Xs = (np.arange(N, dtype=float) - (N / 2) + 1/2) * dx
    elif pp_centering=='odd':
        Xs = (np.arange(N, dtype=float) - (N / 2)) * dx

    du = psf_pixelscale_lamD
    if fp_centering=='odd':
        Us = (np.aarange(npsf, dtype=float) - npsf / 2) * du
    elif fp_centering=='even':
        Us = (np.aarange(npsf, dtype=float) - npsf / 2 + 1/2) * du

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
