import numpy as np
import matplotlib.pyplot as plt

# Simulation Parameters Class
class Simulation:
    def __init__(self, Nx, Ny, dx, dy):
        self.Nx = Nx
        self.Ny = Ny
        self.dx = dx
        self.dy = dy

# Center-Origin Angular Spectrum Method
def angular_spectrum_method_center(simulation, E, z, λ):
    """
    Angular Spectrum Method with center-origin convention.
    """
    n = simulation.Nx
    delkx = 2.0 * np.pi / (simulation.dx * n)
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    k = 2.0 * np.pi / λ

    # Step 1: Apply FFT shift, FFT, and IFFT shift (Center-Origin)
    wf_as = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(E)))

    # Step 2: Generate kx, ky, and kz arrays
    kx, ky = np.meshgrid(kxy, kxy)

    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))

    # Step 3: Transfer function
    tf = np.exp(1j * kz * z)

    # Step 4: Propagate the wavefront and apply inverse FFT
    prop_wf = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(wf_as * tf)))
    return prop_wf



def angular_spectrum_method_center2(simulation, E, z, λ):
    """
    Angular Spectrum Method with center-origin convention.
    """
    n = simulation.Nx
    delkx = 2.0 * np.pi / (simulation.dx * n)
    kxy = np.array([(-(n / 2.0) + i + 0.5) * delkx for i in range(n)])
    k = 2.0 * np.pi / λ

    # Step 1: Apply FFT shift, FFT, and IFFT shift (Center-Origin)
    wf_as = np.fft.fft2(np.fft.fftshift(E))

    # Step 2: Generate kx, ky, and kz arrays
    kx, ky = np.meshgrid(kxy, kxy)
    kz = np.sqrt(np.abs(k**2 - kx**2 - ky**2 + 0j))

    # Step 3: Transfer function
    tf = np.fft.ifftshift(np.exp(1j * kz * z))

    # Step 4: Propagate the wavefront and apply inverse FFT
    prop_wf = np.fft.fftshift(np.fft.ifft2(wf_as * tf))
    return prop_wf



# Main script for comparison
if __name__ == "__main__":
    # Parameters
    Nx, Ny = 1024, 1024  # Grid size
    dx = dy = 1e-6     # Pixel size (1 micron)
    λ = 0.5e-6         # Wavelength (0.5 micron)
    z = 0.01           # Propagation distance (1 cm)

    # Create random input wavefront data
    np.random.seed(42)
    E = np.random.random((Nx, Ny)) + 1j * np.random.random((Nx, Ny))

    # Simulation setup
    sim = Simulation(Nx, Ny, dx, dy)

    # Propagate using both methods
    E_center = angular_spectrum_method_center(sim, E, z, λ)
    E_corner = angular_spectrum_method_center2(sim, E, z, λ)

    # Compute real and imaginary differences
    real_difference = np.real(E_center) - np.real(E_corner)
    imag_difference = np.imag(E_center) - np.imag(E_corner)

    # Plot the results
    plt.figure(figsize=(15, 6))

    # Real components comparison
    plt.subplot(2, 3, 1)
    plt.title("Real Part (Center-Origin)")
    plt.imshow(np.real(E_center), cmap='viridis')
    plt.colorbar()

    plt.subplot(2, 3, 2)
    plt.title("Real Part (Corner-Origin)")
    plt.imshow(np.real(E_corner), cmap='viridis')
    plt.colorbar()

    plt.subplot(2, 3, 3)
    plt.title("Real Part Difference")
    plt.imshow(real_difference, cmap='coolwarm')
    plt.colorbar()

    # Imaginary components comparison
    plt.subplot(2, 3, 4)
    plt.title("Imaginary Part (Center-Origin)")
    plt.imshow(np.imag(E_center), cmap='viridis')
    plt.colorbar()

    plt.subplot(2, 3, 5)
    plt.title("Imaginary Part (Corner-Origin)")
    plt.imshow(np.imag(E_corner), cmap='viridis')
    plt.colorbar()

    plt.subplot(2, 3, 6)
    plt.title("Imaginary Part Difference")
    plt.imshow(imag_difference, cmap='coolwarm')
    plt.colorbar()

    plt.tight_layout()
    plt.show()
