import numpy as np
import time
import math


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


# Parameters
sigma = 10
intensity = 1e5
noise_stddev = 0.01

wavelength = 500e-9          # Wavelength in meters
distance = 1000e-3          # Distance in meters


# Main function to iterate over sizes
def main():
    sizes = [128, 256, 512, 1024, 2048]  # Matrix sizes
    rounds = 25                          # Number of timing rounds

    for mat_rows in sizes:
        mat_size = mat_rows * mat_rows
        pixel_scale = 10e-3 / (mat_rows / 2)  # Adjust pixel scale for each size

        time_arth_mean = 0
        time_geom_mean = 0

        print(f"Processing size: {mat_rows}x{mat_rows}")
        
        for r in range(rounds):
            # Generate Gaussian star pattern
            arr = generate_star_gaussian(mat_rows, sigma, intensity, noise_stddev,True)

            # Time angular spectrum propagation
            start_time = time.time()
            arr = np.fft.fft2(arr)
            end_time = time.time()

            # Calculate elapsed time
            round_time = end_time - start_time
            time_arth_mean += round_time
            time_geom_mean += math.log(round_time)

            # print(f"Round {r+1}: Time = {round_time:.6f} seconds")

        # Finalize mean calculations
        time_arth_mean /= rounds
        time_geom_mean = math.exp(time_geom_mean / rounds)

        # Print results for the current size
        print(f"\nResults for size {mat_rows}x{mat_rows}:")
        print(f"Arithmetic Mean Time: {time_arth_mean:.6f} seconds")
        print(f"Geometric Mean Time: {time_geom_mean:.6f} seconds\n")

if __name__ == "__main__":
    main()
