/*
 * MIT STAR Lab
 * Modified by Subhi in Nov 20, 2024
 */
#include "angular_spectrum_propagation.h"

/**
 * @details
 * The propagation process involves the following steps:
 * 1. **Shift to Frequency Domain**:
 *    - The input wavefront is shifted using a 2D FFT shift to align the zero frequency component 
 *      to the center of the grid.
 *    - A forward 2D FFT is performed to transform the wavefront into the frequency domain.
 * 2. **Apply Transfer Function**:
 *    - A transfer function, representing the angular spectrum propagation kernel, is computed.
 *    - The wavefront in the frequency domain is multiplied element-wise by this transfer function.
 * 3. **Transform Back to Spatial Domain**:
 *    - The result is shifted, and an inverse 2D FFT is applied to transform the data back 
 *      to the spatial domain.
 *    - A final shift aligns the output correctly.
 * 
 * The formula used for the transfer function is:
 *   H(fx, fy) = exp(i * kz * distance),
 * where:
 *   kz = sqrt((2π/λ)^2 - fx^2 - fy^2),
 *   fx and fy are spatial frequencies derived from the grid size and pixel scale.
 * 
 * Python summary:
 * p:    wf_as = xp.fft.ifftshift(xp.fft.fft2(xp.fft.fftshift(wavefront)))
 * p:    kx, ky = xp.meshgrid(kxy,kxy)
 * p:    kz = xp.sqrt(k**2 - kx**2 - ky**2 + 0j)
 * p:    tf = xp.exp(1j*kz*distance.to_value(u.m))
 * p:    prop_wf = xp.fft.fftshift(xp.fft.ifft2(xp.fft.ifftshift(wf_as * tf)))
 * 
 * @note
 * - The scaling of the 2D FFT is included in the transfer function matrix.
 * - Optimiziation is done by pre-shifting the transfer function to remove redundant shifts of the wavefront during multiplication.
 * - A space optimization to make is to save only part of the tranfer function as it's symmetric on all axes
 */
void angular_spectrum_propagation(cmpx_data_t* wavefront, int n, float wavelength, float distance, float pixel_scale)
{
    fftshift2d(wavefront, n, n);
    fft2d(wavefront, n, n, 0, 0);   // invert = 0; scale = 0;

    cmpx_data_t transfer_function[n*n];
    construct_ishifted_transform_function(transfer_function, n, wavelength, distance, pixel_scale);
    elementwise_cmpx_mul(wavefront, transfer_function, n, n);

    fft2d(wavefront, n, n, 1, 0);   // invert = 1; scale = 0; scale included in tf
    fftshift2d(wavefront, n, n);
}

/**
 * @note: this function does everything as doubles then converts to cmpx_data_t
 */
void construct_transform_function(cmpx_data_t* tf, int n, float wavelength, float distance, float pixel_scale)
{
    double wavelength_d  = (double) wavelength;
    double pixel_scale_d = (double) pixel_scale;
    double distance_d    = (double) distance;

    double delkx = 2.0 * M_PI / (pixel_scale_d * n);
    double k = 2.0 * M_PI / wavelength_d;
    double k_2 = k*k;
    double scale = 1/((double)n * (double)n); // 2D FFT normalization

    std::vector<double> kxy(n);
    for (int i = 0; i < n; i++)
        kxy[i] = ((-(double)n / 2.0 + i) + 0.5) * delkx; // Center bins

    int i, j, index;
    double kxy_sum;
    double kz;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            index = i * n + j;
            kxy_sum = kxy[i] * kxy[i] + kxy[j] * kxy[j];
            if (kxy_sum > k_2) {
                tf[index] = 0.0;
                continue;
            }
            kz = sqrt(k_2 - kxy_sum);
            tf[index] = (cmpx_data_t) (scale * std::exp(complex(0.0, kz * distance)));
        }
    }
}



/**
 * @note: this function does everything as doubles then converts to cmpx_data_t
 */
void construct_ishifted_transform_function(cmpx_data_t* tf, int n, float wavelength, float distance, float pixel_scale)
{
    double wavelength_d  = (double) wavelength;
    double pixel_scale_d = (double) pixel_scale;
    double distance_d    = (double) distance;

    double delkx = 2.0 * M_PI / (pixel_scale_d * n);
    double k = 2.0 * M_PI / wavelength_d;
    double k_2 = k*k;
    double scale = 1/((double)n * (double)n); // 2D FFT normalization
    int half_n = n/2;

    std::vector<double> kxy(n);
    for (int i = 0; i < n; i++)
        kxy[i] = ((-(double)n / 2.0 + i) + 0.5) * delkx; // Center bins


    int i, j;
    int top_left_i, top_left_j;
    int top_right_i, top_right_j;
    int bottom_left_i, bottom_left_j;
    int bottom_right_i, bottom_right_j;
    int index;

    double kxy_sum;
    double kz;

    // Split into 2 loops to avoid the possiblity of cache conflicts between upper and lower parts

    // Filling the top half of tf
    // fill top_left with content that is supposed to be in bottom_right
    // fill top_right with content that is supposed to be in bottom_left
    for (i = 0; i < half_n; i++) {
        for (j = 0; j < half_n; j++) {
            top_left_i = i;
            top_left_j = j;

            top_right_i = i;
            top_right_j = j + half_n;

            bottom_left_i = i + half_n;
            bottom_left_j = j;

            bottom_right_i = i + half_n;
            bottom_right_j = j + half_n;

            // fill top_left with content that is supposed to be in bottom_right
            index = top_left_i * n + top_left_j;
            kxy_sum = kxy[bottom_right_i] * kxy[bottom_right_i] + kxy[bottom_right_j] * kxy[bottom_right_j];
            if (kxy_sum > k_2) {
                tf[index] = 0.0;
                continue;
            }
            kz = sqrt(k_2 - kxy_sum);
            tf[index] = (cmpx_data_t) (scale * std::exp(complex(0.0, kz * distance)));

            // fill top_right with content that is supposed to be in bottom_left
            index = top_right_i * n + top_right_j;
            kxy_sum = kxy[bottom_left_i] * kxy[bottom_left_i] + kxy[bottom_left_j] * kxy[bottom_left_j];
            if (kxy_sum > k_2) {
                tf[index] = 0.0;
                continue;
            }
            kz = sqrt(k_2 - kxy_sum);
            tf[index] = (cmpx_data_t) (scale * std::exp(complex(0.0, kz * distance)));
        }
    }


    // Filling the bottom half of tf
    // fill bottom_left with content that is supposed to be in top_right
    // fill bottom_right with content that is supposed to be in top_left
    for (i = 0; i < half_n; i++) {
        for (j = 0; j < half_n; j++) {
            top_left_i = i;
            top_left_j = j;

            top_right_i = i;
            top_right_j = j + half_n;

            bottom_left_i = i + half_n;
            bottom_left_j = j;

            bottom_right_i = i + half_n;
            bottom_right_j = j + half_n;

            // fill bottom_left with content that is supposed to be in top_right
            index = bottom_left_i * n + bottom_left_j;
            kxy_sum = kxy[top_right_i] * kxy[top_right_i] + kxy[top_right_j] * kxy[top_right_j];
            if (kxy_sum > k_2) {
                tf[index] = 0.0;
                continue;
            }
            kz = sqrt(k_2 - kxy_sum);
            tf[index] = (cmpx_data_t) (scale * std::exp(complex(0.0, kz * distance)));

            // fill bottom_right with content that is supposed to be in top_left
            index = bottom_right_i * n + bottom_right_j;
            kxy_sum = kxy[top_left_i] * kxy[top_left_i] + kxy[top_left_j] * kxy[top_left_j];
            if (kxy_sum > k_2) {
                tf[index] = 0.0;
                continue;
            }
            kz = sqrt(k_2 - kxy_sum);
            tf[index] = (cmpx_data_t) (scale * std::exp(complex(0.0, kz * distance)));
        }
    }

}


void elementwise_cmpx_mul(cmpx_data_t* arr1, cmpx_data_t* arr2, int m, int n)
{
    for (int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            arr1[i*m+j] *= arr2[i*m+j];
        }
    }
}

void generate_star_gaussian(cmpx_data_t* arr, int size, double sigma, double intensity, double noise_stddev, bool noise)
{
    // Define the center of the image
    int center_x = size / 2;
    int center_y = size / 2;

    auto add_gaussian = [&](int offset_x, int offset_y, double local_sigma, double local_intensity)
    {
        for (int x=0; x<size; x++){
            for (int y=0; y<size; y++){
                double distance = std::sqrt(std::pow(x - center_x + offset_x, 2) + std::pow(y - center_y + offset_y, 2));
                int idx = x * size + y;
                arr[idx] +=  (cmpx_data_t) (local_intensity * std::exp(-distance * distance / (2 * local_sigma * local_sigma)));
            }
        }
    };

    std::fill(arr, arr + size * size, cmpx_data_t(0.0, 0.0));
    add_gaussian(0, 0, sigma, intensity);
    add_gaussian(size / 8, size / 8, sigma / 2, intensity / 8);


    if (noise){
        // Add random Gaussian noise
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<double> noise_distribution(0.0, noise_stddev);
        for (int x = 0; x < size; ++x) {
            for (int y = 0; y < size; ++y) {
                int idx = x * size + y;
                arr[idx] +=  (cmpx_data_t) noise_distribution(gen);
            }
        }
    }
}
