#include <time.h>
#include <vector>
#include <iostream>
#include "fft.hxx"

/////// Timing Helpers ///////
#define TIMER(label)  timespec label; clock_gettime(CLOCK_MONOTONIC, &label)
#define ELAPSED(b,a)  (double(b.tv_sec - a.tv_sec)*1000000000.0+double(b.tv_nsec-a.tv_nsec))/1000000000.0

using namespace std;

int main() {
    ///////  Testing Data Parameters  ///////
    int DATA_SIZE = 2048;
    float sampling_rate = 0.01;
    std::vector<float> frequencies = {1, 4, 7};
    std::vector<float> amplitudes = {3, 1, 0.5}; 
    std::vector<float> time_range = {0, sampling_rate * DATA_SIZE};
    bool invert = 0;

    std::vector<cmpx_data_t> software_generated_input_data(DATA_SIZE, 0);
    std::vector<cmpx_data_t> software_generated_ouput_data(DATA_SIZE, 0);

    fft_data_generator(software_generated_input_data, software_generated_ouput_data, frequencies, amplitudes, time_range, invert);

    int rounds = 10;
    TIMER(TIME_START);
    for (int i = 0; i < rounds; i++) {
        fft(software_generated_input_data, invert);
    }
    TIMER(TIME_END);

    cout << "Timing FFT on mac with data size " << DATA_SIZE << ", total time for " << rounds << " rounds = " << ELAPSED(TIME_END, TIME_START) << endl;
    cout << "time for one round is = " << ELAPSED(TIME_END, TIME_START) / rounds << endl;

    return 0;
}
