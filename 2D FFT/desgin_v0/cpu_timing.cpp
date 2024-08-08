#include "host/fft.h"
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <sys/syscall.h>
#include <unistd.h>

using namespace std;

// Timing Helpers
#define TIMER(label)  timespec label; clock_gettime(CLOCK_MONOTONIC, &label);
#define ELAPSED(b, a)  (double(b.tv_sec - a.tv_sec) * 1000000000.0 + double(b.tv_nsec - a.tv_nsec)) / 1000000000.0

double fft2d_timer(cmpx_data_t *input_data, cmpx_data_t *output_data, int m, int n, bool invert) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(100.0, 1000.0);

    vector<vector<cmpx_data_t>> output_data_v(m, vector<cmpx_data_t>(n));

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = i * n + j;
            float real = dis(gen);
            float imag = dis(gen);
            input_data[index] = cmpx_data_t(real, imag);
            output_data_v[i][j] = cmpx_data_t(real, imag);
        }
    }

    TIMER(starting_time);
    fft2d(output_data_v, invert);
    TIMER(finishing_time);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = i * n + j;
            output_data[index] = output_data_v[i][j];
        }
    }

    double time_elapsed = ELAPSED(finishing_time, starting_time);

    return time_elapsed;
}

const int n = 1024;
const int m = 1024;
int trials = 10;

cmpx_data_t input_data[n * m];
cmpx_data_t output_data[n * m];

int main() {
    double total_trials_time = 0;
    int total_trials = trials;  // Store the initial number of trials

    while (trials--) {
        total_trials_time += fft2d_timer(input_data, output_data, m, n, false);
    }

    double avg = total_trials_time / total_trials;

    cout << "Timing 2D FFT on matrix with size (" << m << ", " << n << ") for number of trials = " << total_trials << endl;
    cout << "Total time = " << total_trials_time << endl;
    cout << "Average time over one run = " << avg << endl;

    return 0;
}
