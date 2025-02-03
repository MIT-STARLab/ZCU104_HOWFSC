
#include "angular_spectrum_propagation.h"
#include <sys/syscall.h>
#include <time.h>
#include <complex>
#include <cmath>

using namespace std;

#define TIMER(label)  timespec label; clock_gettime(CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0


int main() {
    // Star parameters
    float sigma = 10;
    float intensity = 1e5;
    float noise_stddev = 0.01;

    // Telescope parameters
    float wavelength = 500e-9;
    float distance = 1000e-3;

    // Matrix sizes to test
    const int sizes[] = {128, 256, 512, 1024, 2048};
    const int num_sizes = sizeof(sizes) / sizeof(sizes[0]);

    // Loop over different matrix sizes
    for (int s = 0; s < num_sizes; s++) {
        int mat_rows = sizes[s];
        int mat_size = mat_rows * mat_rows;
        float pixel_scale = 10e-3 / (mat_rows / 2); // Adjust pixel scale for each size

        // Allocate memory for the complex array
        complex<float> *arr = (complex<float> *)malloc(mat_size * sizeof(complex<float>));
        if (!arr) {
            cerr << "Failed to allocate memory for size " << mat_rows << "x" << mat_rows << "!" << endl;
            continue;
        }

        // Timing variables
        double time_arth_mean = 0;
        double time_geom_mean = 0;

        int rounds = 25; // Number of rounds for averaging
        for (int r = 0; r < rounds; r++) {
            // Generate Gaussian star pattern
            generate_star_gaussian(arr, mat_rows, sigma, intensity, noise_stddev, true);

            // Time the angular spectrum propagation
            TIMER(starting_point);
            angular_spectrum_propagation(arr, mat_rows, wavelength, distance, pixel_scale);
            TIMER(ending_point);

            // Calculate elapsed time
            double round_time = ELAPSED(ending_point, starting_point);
            time_arth_mean += round_time;            // Add to arithmetic mean
            time_geom_mean += log(round_time);       // Add to geometric mean (logarithmic)

            // cout << "Size " << mat_rows << "x" << mat_rows 
            //      << ", Round " << r + 1 << ": Time = " << round_time << " seconds" << endl;
        }

        // Finalize timing calculations
        time_arth_mean /= rounds; // Arithmetic mean
        time_geom_mean = exp(time_geom_mean / rounds); // Geometric mean

        // Output results for the current size
        cout << "\nResults for Size " << mat_rows << "x" << mat_rows << ":\n";
        cout << "Arithmetic Mean Time: " << time_arth_mean << " seconds\n";
        cout << "Geometric Mean Time: " << time_geom_mean << " seconds\n\n";

        // Free allocated memory
        free(arr);
    }

    return EXIT_SUCCESS;
}