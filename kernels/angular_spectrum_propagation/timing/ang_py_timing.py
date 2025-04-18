import numpy as np
import time
import math
from angular_spectrum_testing import *

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
            arr = angular_spectrum_progagation_float(arr, mat_rows, wavelength, distance, pixel_scale)
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
